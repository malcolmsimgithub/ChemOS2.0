import click
import coloredlogs

from . import admin as _admin

@click.group()
def cli():
    pass


@cli.group(help='Admin operations on the database')
def admin():
    pass    
    


@cli.command()
def get():
    # Options -> no_filter, csv (output a csv)
    # Choose table       <-+ 
    # Add another table? --+ (show options)
    # Add filter                       <------------------+
    #   Filter on table -> choose table if more than one  |
    #   Choose field                                      |
    #   Choose operator                                   |
    #   Choose value                                      |
    # Add another filter ---------------------------------+
    # Return results

    pass


@cli.command()
def add():
    # Choose type -> add new
    # Fill data                          
    # Fill relationships
    #   -> relationship does not exists -> add new
    #   -> relationship exists -> get if not too big, get with filter otherwise
    pass


@cli.command()
def update(type, uuid):
    # give uuid and type, edit fields
    # fetch data, print data
    # choose column to edit <-+
    # edit   -----------------+
    # show and confirm
    # update
    pass


@cli.command()
def delete(type, uuid):
    # give uuid and type, delete field
    # fetch, show, confirm, delete
    pass


@admin.command(help='Create a database')
@click.option('-u', '--user', prompt=True)
@click.option('-p', '--password', prompt=True, hide_input=True)
@click.option('-h', '--hostname', prompt=True)
@click.option('-d', '--db_name', prompt=True)
@click.option('-s', '--struct_file', prompt=True)
@click.option('-e', '--es_file', prompt=True)
@click.option('-c', '--config', help="Alembic configuration file")
@click.option('--admin_pass', prompt=True, hide_input=True, confirmation_prompt=True)
@click.option('--user_pass', prompt=True, hide_input=True, confirmation_prompt=True)
def create_db(user, password, hostname, db_name, struct_file, es_file, config, admin_pass, user_pass):
    _admin.create_db(user, password, hostname, db_name, struct_file, es_file, admin_pass, user_pass, 'alembic.ini')


@admin.command(help='Remove a database')
@click.option('-u', '--user', prompt=True)
@click.option('-p', '--password', prompt=True, hide_input=True)
@click.option('-h', '--hostname', prompt=True)
@click.option('-d', '--db_name', prompt=True)
@click.confirmation_option(prompt='Are you sure you want to drop the db?')
def drop_db(user, password, hostname, db_name):
    _admin.drop_db(user, password, hostname, db_name)


@admin.command(help='Backup a database (save the alembic version)')
@click.option('-u', '--user', prompt=True)
@click.option('-p', '--password', prompt=True, hide_input=True)
@click.option('-h', '--hostname', prompt=True)
@click.option('-d', '--db_name', prompt=True)
@click.option('-c', '--compression')
@click.option('-o', '--output_file', type=click.File('wb'), prompt=True)
def backup_db(user, password, hostname, db_name, compression, output_file):
    _admin.backup_db(user, password, hostname, db_name, output_file, compression)


@admin.command(help='Empty a database')
@click.option('-u', '--user', prompt=True)
@click.option('-p', '--password', prompt=True, hide_input=True)
@click.option('-h', '--hostname', prompt=True)
@click.option('-d', '--db_name', prompt=True)
@click.confirmation_option(prompt="Are you sure you want to empty the db?")
def empty_db(user, password, hostname, db_name):
    _admin.empty_db(user, password, hostname, db_name)


@admin.command(help='Restore a database with a dump file')
@click.option('-u', '--user', prompt=True)
@click.option('-p', '--password', prompt=True, hide_input=True)
@click.option('-h', '--hostname', prompt=True)
@click.option('-d', '--db_name', prompt=True)
@click.option('-f', '--dump_file', type=click.File('r'))
@click.confirmation_option(prompt="Are you sure you want to do that?")
def restore_db(user, password, hostname, db_name, dump_file):
    _admin.restore_db(user, password, hostname, db_name, 'alembic.ini', dump_file)



if __name__ == '__main__':
    coloredlogs.install('INFO')
    cli()
