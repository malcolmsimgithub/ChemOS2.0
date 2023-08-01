from .client import MDBClient, MDBClientWithSSH

import toml

__version__ = '0.2'
__author__ = 'Theophile Gaudin'

def load_client_from_config(config_file):
    """
    Loads a configuration file and initialize a client.

    :param config_file: path of the configuration file (string):
    
    :returns: MDBClient

    """
    with open(config_file, 'r') as f:
        config = toml.load(f)
    
    kwargs = config['database']
    kwargs['use_tqdm'] = config['use_tqdm']

    if config['use_ssh']:
        ssh = config['ssh']
        kwargs['ssh_username'] = ssh['username']
        kwargs['ssh_hostname'] = ssh['hostname']
        kwargs['ssh_keyfile'] = ssh['keyfile']
        return MDBClientWithSSH(**kwargs)
    return MDBClient(**kwargs)
