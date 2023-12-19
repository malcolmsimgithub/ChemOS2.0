import sqlalchemy as sa
from alembic.config import Config
from alembic import command

from datetime import datetime
import os
import socket
import logging
import subprocess
import shlex

from .. import MDBClient


def create_connection(user, password, hostname, database):
    return sa.create_engine(f'postgresql://{user}:{password}@{hostname}/{database}').connect()


def create_user(conn, user, password, logger=None):
    logger = logger or logging.getLogger(__name__)
    logger.info(f"Adding {user} to the database")
    conn.execute(f'CREATE ROLE {user} WITH LOGIN PASSWORD \'{password}\';')


def transfer_ownership(conn, database, new_owner, logger=None):
    logger = logger or logging.getLogger(__name__)
    logger.info(f'Transfering ownership to {new_owner}')
    conn.execute(f'''
    DO $$DECLARE r record; 
    BEGIN 
        FOR r in SELECT table_name FROM information_schema.tables WHERE table_schema = \'public\' 
        LOOP 
            EXECUTE \'ALTER TABLE public.\'|| quote_ident(r.table_name) ||\' OWNER TO {new_owner};\'; 
        END LOOP; 
    END$$;
    ''') 

    conn.execute(f'ALTER TABLE sourcing.eventstore OWNER TO {new_owner};')
    conn.execute(f'ALTER SCHEMA public OWNER TO {new_owner};')
    conn.execute(f'ALTER SCHEMA sourcing OWNER TO {new_owner};')
    conn.execute(f'ALTER DATABASE {database} OWNER TO {new_owner};')
    conn.execute(f'ALTER FUNCTION public.create_synthesis_hid OWNER TO {new_owner};') 
    conn.execute(f'ALTER FUNCTION sourcing.on_event OWNER TO {new_owner};')


def grant_access_right(conn, user, logger=None):
    logger = logger or logging.getLogger(__name__)
    logger.info(f'Granting usage access to {user}')
    conn.execute(f'GRANT USAGE ON SCHEMA public TO {user};')
    conn.execute(f'GRANT USAGE ON SCHEMA sourcing TO {user};')
    conn.execute(f'GRANT EXECUTE ON FUNCTION public.create_synthesis_hid TO {user};')
    conn.execute(f'GRANT EXECUTE ON FUNCTION sourcing.on_event TO {user};')
    conn.execute(f'GRANT USAGE ON ALL SEQUENCES IN SCHEMA public TO {user};')
    conn.execute(f'GRANT USAGE ON ALL SEQUENCES IN SCHEMA sourcing TO {user};')
    conn.execute(f'GRANT SELECT, INSERT ON TABLE sourcing.eventstore TO {user};')
    conn.execute(f'GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO {user};') 


def create_db(user, password, hostname, new_db_name, structure_file,
        event_sourcing_file, admin_password, user_password, alembic_conf,
        alembic_version='head'):
    logger = logging.getLogger(__name__)
    
    logger.info(f"Creating database {new_db_name}.")
    conn = create_connection(user, password, hostname, 'postgres')
    conn.execution_options(isolation_level="AUTOCOMMIT").execute(f"CREATE DATABASE {new_db_name};")
    conn.close()

    conn = create_connection(user, password, hostname, new_db_name)
    
    structure = sa.text(open(structure_file, 'r').read())
    event_sourcing = sa.text(open(event_sourcing_file, 'r').read())

    logger.info("Loading database schema and event sourcing")
    conn.execute('CREATE EXTENSION "uuid-ossp";') 
    conn.execute(structure)
    conn.execute(event_sourcing)

    admin_username = new_db_name +'_admin'
    create_user(conn, new_db_name, user_password, logger)
    create_user(conn, admin_username, admin_password, logger)

    # Transferring ownership to admin 
    transfer_ownership(conn, new_db_name, admin_username, logger)
    
    logger.info("Performing alambic updates")
    conf = get_alembic_conf(
        f'postgresql://{new_db_name}_admin:{admin_password}@{hostname}/{new_db_name}',
        alembic_conf
    )
    alembic_upgrade(conf, alembic_version)

    # Granting rights
    grant_access_right(conn, new_db_name, logger)
    conn.close()
    logger.info("Done!")


def alembic_upgrade(conf, version, conn=None):
    if conn is not None:
        conf.attributes['connection'] = conn
    command.upgrade(conf, version)


def drop_db(user, password, hostname, db_name):
    conn = create_connection(user, password, hostname, 'postgres')
    logger = logging.getLogger(__name__)
    logger.info(f'Dropping database {db_name}')
    conn.execution_options(isolation_level="AUTOCOMMIT").execute(f'DROP DATABASE {db_name};')
    conn.execute(f'DROP ROLE {db_name};')
    conn.execute(f'DROP ROLE {db_name}_admin;')
    conn.close()
    logger.info("Done!")


def get_alembic_conf(sql_url, alembic_conf):
    conf = Config(alembic_conf)
    conf.set_main_option('sqlalchemy.url', sql_url)
    return conf


def backup_db(user, password, hostname, db_name, output_file, compression, logger=None):
    logger = logger or logging.getLogger(__name__)
    out = subprocess.check_output(['which', 'pg_dump'])
    if out.decode().endswith('not found'):
        logger.warn("pg_dump has not been found! Please install postgresql-client.")
        return

    client = MDBClient(hostname, user, password, db_name)
    alembic_version = client.session.query(client.models.alembic_version).one().version_num

    output_file.write(f'-- alembic_version={alembic_version}\n'.encode())
    output_file.write((f'-- backup of {hostname}/{db_name} by '
                       f' {user} @ {socket.gethostname()} on '
                       f'{datetime.now()}\n').encode())
    
    command = (f'pg_dump --data-only -Z {compression} --no-owner '
               f' --host {hostname} -d {db_name} -U {user} '
               f' -t sourcing.eventstore')

    logger.info(f'Backing up {db_name} @ {hostname} in {output_file.name}')
    os.environ["PGPASSWORD"] = password
    subprocess.Popen(shlex.split(command), stdout=output_file)


def empty_db(user, password, hostname, db_name):
    """
    This empty the eventstore, but does not reset the indexing.

    """
    client = MDBClient(hostname, user, password, db_name)
    client.rollback(datetime(1950, 1, 1))
    client.session.query(client.models.eventstore).delete()
    client.session.commit()


def restore_db(user, password, hostname, db_name, alembic_conf, dump_file):
    out = subprocess.check_output(['which', 'pg_dump'])
    if out.decode().endswith('not found'):
        print("pg_restore has not been found! Please install postgresql-client.")
        return

    empty_db(user, password, hostname, db_name)

    alembic_version = dump_file.readline().split('=')[-1]

    conf = get_alembic_conf(f'postgresql://{user}:{password}@{hostname}/{db_name}',
                            alembic_conf)
    alembic_upgrade(conf, alembic_version)

    os.environ['PGPASSWORD'] = password
    command = (f'pg_restore -x --data-only --if-exists --no-owner -t '
               f' sourcing.eventstore -h {hostname} -U {user} -d {db_bame}'
               f' {dump_file.name} ') 
