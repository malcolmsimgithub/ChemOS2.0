import pandas as pd
import numpy as np
from tqdm import tqdm
import logging
import sqlalchemy
from sqlalchemy import and_, text
from threading import Thread
import json

from . import ssh
from . import database as db
from . import mapper


class MDBClient(mapper.SchemaMapper):
    """
    MDBClient is a high-level client for the madness database. There is an
    emphasis on user-friendliness and a thight integration with pandas.

    :param hostname: database hostname (string)
    :param username: database username (string)
    :param password: database password (string)
    :param database: name of the database to access (string)
    :param port: port to use to connect to the database (integer, default: 5432)
    :param use_tqdm: how a tqdm progressbar if set to True. (boolean, optional)

    """
    def __init__(self, hostname, username, password, database, port=5432, use_tqdm=True, logger=None):
        self.sql_url = f"postgresql://{username}:{password}@{hostname}:{port}/{database}"
        self.logger = logger or logging.getLogger(__name__)
        Session, engine, models = db.init_db(self.sql_url)
        self.session = Session()
        self.engine = engine
        self.models = models

        self.dao =  DataAccessObject(self.session, self.models)
        self.rollback = self.dao.rollback

        self.use_tqdm = use_tqdm

    def __delete__(self):
        self.session.close()

    def get(self, table_name, return_df=True, filters=None, limit=None, offset=None, order_by=None):
        """
        This method allows to read data stored in the database.

        :param table_name: (str or List(str)) name(s) of the table(s) to access
        :param return_df: (bool) returns a pandas.DataFrame if True
        :param filters: (list) list of fitlers 
        :param limit: (int or None) limits how many rows should be loaded.
        :param offset: (int or None) indicates how many rows to skip.
        :param order_by: (asc or desc clause) order the query. The default is
        by timestamp or created_on column.
        
        :returns: List of queried object of pandas.DataFrame containing the data

        Typical usage:        
       
        .. code-block:: python

            client = MDBCLient('localhost', 'postgres', '', 'madness')
            df = client.get('fragment')
            df.head()

        It is also possible to acces data from many tables at the same time, as
        long as they are related. For instance, in the here-below example, we
        are joining both table `fragment` and `molecule` and filter for a
        specific smiles.
        
        .. code-block:: python

            client = MDBCLient('localhost', 'postgres', '', 'madness')
            df = client.get(['fragment', 'molecule'],
                            filters=[client.models.Molecule.smiles == 'C1=CC=CC=C1'])
            df.head()
        
        """
        data = self.dao.get(table_name, filters, limit, offset, order_by)
        if not return_df:
            return data
        df = pd.DataFrame.from_records([d.__dict__ for d in data])
        if '_sa_instance_state' in df.columns:
            del df['_sa_instance_state']
            
        df = df.replace({np.nan: None})
        return df
    
    def get_id(self, table_name, **kwargs):
        """
        This method helps finding the id of a specific item in the database.
        it will be usefule when adding data with relationships.

        :param table_name (str): name of the table
        :param kwargs: keyword arguments which will be used to build filters. Those 
        keywords have to correspond to existing column in the table.

        :returns: The id as a string. If the objects can't be find, a
        exception is raised.

        .. code-block:: python
            
            client = MDBCLient('localhost', 'postgres', '', 'madness')
            molecule_id = client.get_id('molecule', smiles='C1=CC=CC=C1')

        """
        model = getattr(self.models, table_name)
        filters = []
        for k, v in kwargs.items():
            field = getattr(model, k)
            if field is None:
                raise ValueError(f'{table_name} has no column {k}!')
            filters.append(field == v)

        q = self.get(table_name, filters=filters, return_df=False, limit=1).one_or_none()
        if q is None:
            raise ValueError(f'No row with the specified filters where found.  Filters: {filters}')
        return getattr(q, f'{table_name}_id')

    def add(self, table_name, data):
        """
        This method can add data to any existing table in the database. This
        method will returns the events created in the eventstore which contains
        some information about the object such as its id.

        It can be used to perform "batch" operation, although it will fail if a
        constraint is violated.

        :param table_name: name of the table (string)
        :param data: data to add (list, pandas.DataFrame or dict)

        Example:

        .. code-block:: python
            
            client = MDBCLient('localhost', 'postgres', '', 'madness')
            event = client.add('molecule', {'smiles': 'C1=CC=CC=C1'})
            print(event.uuid)
        
        """
        if isinstance(data, dict):
            data = [data]
        if isinstance(data, list):
            data = pd.DataFrame.from_records(data)
        if not isinstance(data, pd.DataFrame):
            raise NotImplementedError("Incorrect data types")

        data = data.replace({np.nan: None})
        
        iter = data.iterrows()
        if self.use_tqdm:
            iter = tqdm(iter, total=len(data))

        events = []
        for i, record in enumerate(iter):
            rec = record[1].replace({pd.np.nan: None})
            events.append(self.dao.add(table_name, rec.to_dict()))
            if i % 500 == 0:
                self.session.commit()
        self.session.commit()
        return events

    def update(self, table_name, data, id=None):
        """
        This method serves to update specifics rows of a table.

        :param table_name (str): name of the table
        :param data (pandas.DataFrame, list or dict): data to update
        :param id (str or list): if data does not contains the id, it must be specified.

        Example:

        .. code-block:: python
            
            client = MDBClient('localhost', 'postgres', '', 'madness')
            df = client.get('molecule')
            df.at[0, 'smiles'] = 'C#N'
            client.update('molecule', df)
            df = client.get('molecule')
            df.head()
        """
        if isinstance(data, dict):
            data = [data]
        if isinstance(data, list):
            exclude = [x for x in data[0].keys() \
                          if x == f'{table_name}_id' or \
                             x == 'created_on' or \
                             x == 'updated_on']
            data = pd.DataFrame.from_records(data, exclude=exclude)
        if not isinstance(data, pd.DataFrame):
            raise TypeError((f"Invalid data type! Got {type(data)} "
                              "but expected either dict, list or "
                              "pandas.DataFrame"))

        data = data.replace({pd.np.nan: None})

        if not id:
            id = data[f'{table_name}_id'].tolist()
            del data[f'{table_name}_id']
        
        if not isinstance(id, list):
            id = [id]
        
        if 'updated_on' in data:
            del data['updated_on']

        if 'created_on' in data:
            del data['created_on']
        
        iter = data.iterrows()
        if self.use_tqdm:
            iter = tqdm(iter)

        events = []
        for i, (record, id) in enumerate(zip(iter, id)):
            rec = record[1].replace({pd.np.nan: None})
            events.append(self.dao.update(table_name, rec.to_dict(), id))
            if i % 500 == 0:
                self.session.commit()
        self.session.commit()
        return events

    def delete(self, table_name, id):
        """
        This method deletes data from the database.

        :param table_name (str): name of the table
        :param id (str, List(str)): id to remove

        .. code-block:: python
            
            client = MDBClient('localhost', 'postgres', '', 'madness')
            df = client.get('molecule')
            client.delete('molecule', df['id'].tolist())
            

        """
        if isinstance(id, str):
            id = [id]

        events = []
        for i, id in enumerate(id):
            events.append(self.dao.delete(table_name, id))
            if i % 500 == 0:
                self.session.commit()
        self.session.commit()
        return events

    def rollback(self, before):
        """
        Performs a rollback on the database. For now, it only supports to be
        restored to a prior date.

        :param before (datetime.datetime): prior date.
        
        Example:

        .. code-block:: python
        
            from datetime import datetime
            client = MDBClient('localhost', 'postgres', '', 'madness')
            client.rollback(datetime(1980, 12, 25)
        """

        self.dao.rollback(before)


class MDBClientWithSSH(MDBClient):
    """
    This version of the client creates an SSH tunnel before connecting to the
    database. This is useful if you are not in the same network as the
    database.

    """
    def __init__(self, hostname, username, password, database,
            ssh_username, ssh_hostname, ssh_keyfile, use_tqdm=True):
        transport = ssh.connect(ssh_username, ssh_hostname, ssh_keyfile)
        port = ssh.get_free_tcp_port()
        forward_server = ssh.tunnel_factory(port, hostname, 5432, transport)
        server_thread = Thread(target=forward_server.serve_forever)
        server_thread.daemon = True
        server_thread.start()
        self.ssh_server = forward_server

        super().__init__('localhost', username, password, database, port)

    def __delete__(self):
        super().__delete__()
        self.ssh_server.shutdown()


class DataAccessObject:
    """
    Defines the low level interactions between MDBClient and the database. It is
    not meant to be used on its own.

    :session: (sqlalchemy.orm.session.Session)
    :models: (sqlalchemy.orm.mapper)
    """
    def __init__(self, session, models):
        self.session = session
        self.models = models

    def get(self, table_name, filters=None, limit=None, offset=None,
            order_by=None):
        """
        Reads data from a table

        :param table_name (str):
        :param filters (list):
        :param limit (int or None):
        :param offset (int or None):
        :param order_by (asc or desc clause):
        """

        if not isinstance(table_name, list):
            table_name = [table_name]

        t = table_name.pop(0)
        model = getattr(self.models, t, None)
        assert model is not None, f'{t} does not correspond to any table!'

        order_by_col = getattr(model, 'updated_on', None)
        if order_by_col is None:
            order_by_col = getattr(model, 'timestamp', None)

        if order_by is None:
            order_by = order_by_col.desc()
        query = self.session.query(model)

        for t in table_name:
            model = getattr(self.models, t, None)
            assert model is not None, f'{t} does not correspond to any table!'
            query = query.join(model)

        if filters is not None:
            filter = True
            for f in filters:
                filter = and_(filter, f)
            query = query.filter(filter)
        query = query.order_by(order_by)

        if limit is None:
            return query.all()

        if offset is not None:
            query = query.offset(offset)
        query = query.limit(limit)

        return query

    def add(self, table_name, data):
        """
        Adds data to the database. This corresponds to a 'create' event on the
        eventstore with a 'type' of table_name

        :param table_name:
        :param data:
        """
        params = {'event': 'create',
                  'type': table_name,
                  'data': data}
        event = self.models.eventstore(**params)
        self.session.add(event)
        return event

    def delete(self, table_name, id):
        """
        Deletes data from the databse. This corresponds to a 'delete' event on
        the eventstore with a 'type' of table_name

        :param table_name:
        :param id:
        """
        params = {'event': 'delete',
                  'type': table_name,
                  'uuid': id}
        event = self.models.eventstore(**params)
        self.session.add(event)
        return event

    def update(self, table_name, data, id):
        """
        Updates data from the database. This corresponds to an 'update' event on
        the eventstore with a 'type' of table_name

        :param table_name:
        :param data:
        :param id:
        """
        params = {'event': 'update',
                  'type': table_name,
                  'data': data,
                  'uuid': id}
        event = self.models.eventstore(**params)
        self.session.add(event)
        return event

    def rollback(self, before):
        """
        Rollback the database. This corresponds to a 'rollback' event on the
        eventstore. For now, until rollback to a specific date is available.

        :param before: (datetime.dateime)
        """

        params = {'event': 'rollback',
                  'data': {'before': before.isoformat()}}
        event = self.models.eventstore(**params)
        self.session.add(event)
        return event

    def commit_or_fetch_event(self, table_name, data):
        """
        Tries to commit the current transation and if it fails,
        fetches the corresponding event in the eventstore table.
        """
        try:
            self.session.commit()
        except sqlalchemy.exc.IntegrityError as ex:
            if 'UniqueViolation' in ex._message():
                self.session.rollback()
                filters = ''
                for k, v in data.items():
                    if isinstance(v, dict):
                        filters += f"data->'{k}' = '{json.dumps(v)}'::jsonb and "
                        continue
                    filters += f"data->>'{k}' = '{v}' and "

                # Single query to retrieve the event
                sql = text(f"""
                    select uuid,
                           type,
                           jsonb_agg(data)->>-1 as data,
                           string_agg("event", ',') as event,
                           max(timestamp) as timestamp
                      from sourcing.eventstore
                      join {table_name}
                        on eventstore.uuid = {table_name}_id
                     where {filters[:-4]}
                  group by uuid, type
                    having not right(string_agg("event", ','), 6) = 'delete'
                  order by min(eventstore.id) asc;
                """)
                out = self.session.execute(sql).fetchall()
                # It's a UniqueViolation... there should be only one
                if len(out) != 1:
                    raise ex

                # Returning uuid as str
                out[0]._processors[0] = lambda x: str(x)
                return out[0]
            else:
                self.session.rollback()
                raise ex
