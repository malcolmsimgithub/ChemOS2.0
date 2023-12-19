from mdb import MDBClient
from mdb import MDBClientWithSSH
#from . import mdb 
import numpy as np
import json
import matplotlib.pyplot as plt


class MDB_client_TE():
    """
    MDB_client_TE is a high-level client of the madness database for the Transient PL measurement.
    """ 

    def __init__(self, config_file, admin = False):

        self.config_file = config_file 
        self._open_config()

        self.client = MDBClient(hostname=self.config['database']['hostname'], 
                                username=self.config['database']['username'], 
                                password=self.config['database']['password'],
                                database=self.config['database']['database'])

        if not admin:
            self._parse_config()

    def _open_config(self):
        with open(self.config_file) as content:
            self.config = json.loads(content.read())

    def _parse_config(self):
        # parse configuration

        self.lab_id = self.get_id('lab', name = self.config['lab']['name'])
        self.machine_id = self.get_id('experiment_machine', name = self.config['experiment_machine']['name'])
        self.exp_type_id = self.get_id('experiment_type', name = self.config['experiment_type']['name'])
        self.unit_ids = {key : self.get_id('data_unit', name = val) for key, val in self.config['units'].items()}


    def _apply_config(self):

        if self.config['experiment_machine']['metadata'] == "None": 
            self.config['experiment_machine']['metadata'] = None
        #add lab and experiment_type
        self.client.add_lab(self.config['lab']['name'], self.config['lab']['short_name'])
        self.client.add_experiment_type(self.config['experiment_type']['name'])

        lab_id = self.client.get_id('lab', name = self.config['lab']['name'])
        experiment_type_id = self.client.get_id('experiment_type', name = self.config['experiment_type']['name'])
        #add machine
        if len(self.get('experiment_machine', \
               filters = [self.client.models.experiment_machine.name == self.config['experiment_machine']['name']])) == 0:
            self.client.add_experiment_machine(self.config['experiment_machine']['name'], self.config['experiment_machine']['make'],\
                self.config['experiment_machine']['model'], self.config['experiment_machine']['metadata'] , experiment_type_id, lab_id)
        else : 
            print('experiment_machine %s already exist' %self.config['experiment_machine']['name'])
        #add units
        for val in self.config['units'].values():
            self.client.add_data_unit(val)

        self._parse_config()
        

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
        return self.client.get_id(table_name, **kwargs)


    def get(self, table_name, filters=None, limit=None, offset=None, order_by=None):
        """
        Reads data from a table

        :param table_name (str): 
        :param filters (list):
        :param limit (int or None):
        :param offset (int or None):
        :param order_by (asc or desc clause):

        :returns: List of queried object of pandas.DataFrame containing the data

        .. code-block:: python

            df = client.get(['fragment', 'molecule'],
                            filters=[client.models.Molecule.smiles == 'C1=CC=CC=C1'])
            df.head()
        """
        return self.client.get(table_name, filters=filters, limit=limit, offset=offset, order_by=order_by)


    def get_local_experiments(self):
        """
        get experiments done by this machine.

        :param filters (list):

        :returns: List of queried object of pandas.DataFrame containing the data
        """
        return self.get('experiment', \
               filters = [self.client.models.experiment.experiment_machine_id == self.machine_id])


    def get_experiment(self, experiment_id = None, synthesis_id = None, molecule_id = None):
        """
        get experiments specified by experiment_id or sysnthesis_is or molecule_id.

        :param experiment_id (str or None):
        :param synthesis_id (str or None):
        :param molecule_id (str or None):

        :returns: List of queried object of pandas.DataFrame containing the data
        """
        if sum([experiment_id is not None, synthesis_id is not None, molecule_id is not None]) != 1:
            raise Exception(' experiment_id xor synthesis_id xor molecule_id should be specified.')

        if experiment_id:
            return self.get('experiment', \
                filters = [self.client.models.experiment.experiment_id == experiment_id])
        if synthesis_id:
            return self.get('experiment', \
                filters = [self.client.models.experiment.synthesis_id == synthesis_id])
        if molecule_id:
            return self.get('experiment', \
                filters = [self.client.models.experiment.molecule_id == molecule_id])



    def get_unit_name(self, data_unit_id):
        """
        get unit name.

        :param data_unit_id(uuid):

        :returns: The name of the data_unit as a string. 
        """
        return self.get('data_unit', \
               filters = [self.client.models.data_unit.data_unit_id == data_unit_id])['name'][0]


    def delete(self, table_name, id):
        """
        Deletes data from the databse. This corresponds to a 'delete' event on
        the eventstore with a 'type' of table_name

        :param table_name(string):
        :param id(uuid):
        """
        return self.client.delete(table_name, id)


    def delete_experiment(self, experiment_id):
        """
        delete experiment with the associated xy/xyz data and child experiments.

        :param experiment_id (uuid):
        """
        def _delete_data(exp_id):
            xy_data = self.get('xy_data', filters = {self.client.models.xy_data.experiment_id ==exp_id})
            xyz_data = self.get('xyz_data', filters = {self.client.models.xyz_data.experiment_id ==exp_id})
            
            if len(xy_data) != 0:
                for xy_id in xy_data['xy_data_id']:
                    self.delete('xy_data', xy_id)
            if len(xyz_data) != 0: 
                for xyz_id in xyz_data['xyz_data_id']:
                    self.delete('xyz_data', xyz_id)

        child_experiments = self.get('experiment', filters = {self.client.models.experiment.parent_experiment_id ==experiment_id})

        if len(child_experiments) != 0: 
            for child_experiment_id in child_experiments['experiment_id']:
                _delete_data(child_experiment_id)
                self.delete('experiment', child_experiment_id)

        _delete_data(experiment_id)

        self.delete('experiment', experiment_id)


    def update(self, table_name, data, id):
        """
        Updates data from the database. This corresponds to an 'update' event on
        the eventstore with a 'type' of table_name

        :param table_name (str): name of the table
        :param data (pandas.DataFrame, list or dict): data to update
        :param id (str or list): id of the data to update

        Example:

        .. code-block:: python
            client.update('experiment', {'notes': 'test experiment'}, experiment_id)
        """
        return self.client.update(table_name, data, id)


    def add_metadata(self, table_name, item_id, metadata):
        """
        add metadata to the item

        :param table_name(str): 
        :param item_id(uuid): 
        :param metadata(dict):

        returns: the events created in the eventstore
        """
        return self.update(table_name, {'metadata' : metadata }, id = item_id)


    def add_notes(self, table_name, item_id, notes):
        """
        add notes to the item

        :param table_name(str): 
        :param item_id(uuid):
        :param notes(text):

        returns: the events created in the eventstore
        """
        return self.update(table_name, {'notes' : notes }, id = item_id)


    def add_xy_data(self, experiment_id, data_name, data, x_units_id, y_units_id, metadata = None, notes = None):
        """
        add xy to the xy_data table.

        :param experiment_id(uuid): experiment_id to which the xy_data to add belongs.
        :param data_name(str): name of the xy data 
        :param data (list or numpy array) : xy data to add([[x_list][y_list]])
        :param x_units_id(uuid): id of the data_units for x
        :param y_units_id(uuid): id of the data_units for y
        :param metadata(dict): metadata of the xy_data
        :param notes(text): notes of the xy_data

        returns: the events created in the eventstore
        """
        for i, d in enumerate(data):
            if type(d) != list:
                data[i] = d.tolist()

        xy_data = self.client.add_xy_data_experiment(experiment_id = experiment_id,
                                                    name = data_name,
                                                    x = data[0],
                                                    y = data[1],
                                                    x_units_id = x_units_id,
                                                    y_units_id = y_units_id)
        if metadata:
            self.add_metadata('xy_data', xy_data.uuid, metadata)
        if notes:
            self.add_notes('xy_data', xy_data.uuid, notes)

        return xy_data


    def add_xyz_data(self, experiment, data_name, data, x_units_id, y_units_id, z_unit_id, metadata = None, notes = None):
        """
        add xyz to the xyz_data table.

        :param experiment_id(uuid): experiment_id to which the xyz_data to add belongs.
        :param data_name(str): name of the xyz data 
        :param data (dict) : xyz data to add ({'x' : [], 'y' : [], 'z': [[][]...,]})
        :param x_units_id(uuid): id of the data_units for x
        :param y_units_id(uuid): id of the data_units for y
        :param z_units_id(uuid): id of the data_units for z
        :param metadata(dict): metadata of the xy_data
        :param notes(text): notes of the xy_data

        returns: the events created in the eventstore
        """
        for key, val in data.items():
            if type(val) != list:
               data[key] = val.tolist()

        xyz_data = self.client.add_xyz_data_experiment(experiment_id = experiment.uuid,
                                           name = data_name,
                                           x = data['x'],
                                           y = data['y'],
                                           z = data['z'],
                                           x_units_id = x_units_id,
                                           y_units_id = y_units_id,
                                           z_units_id = z_units_id)
        if metadata:
            self.add_metadata('xyz_data', xyz_data.uuid, metadata)
        if notes:
            self.add_notes('xyz_data', xyz_data.uuid, notes)

        return xyz_data


    def add_experiment(self, metadata, raw_data_path, synthesis_id = None, molecule_id = None, parent_experiment_id = None, notes = None):
        """
        add new experiment to the experiment table.

        :param metadata(dict): metadata of the xy_data
        :param raw_data_path(str): data path of the raw experiment data
        :param synthesis_id(uuid or None): id of the synthesis which the experiment to add associated with.
        :param molecule_id(uuid or None): id of the molecule targeted in the experiment.
        :param parants_experiment_id(uuid or None): id of the experiment which the experiment to add associated with.
        :param notes(text): notes of the xy_data.

        returns: the events created in the eventstore
        """
        if (synthesis_id is None and molecule_id is None) or \
           (synthesis_id is not None and molecule_id is not None):
            raise ValueError('Either synthesis_id xor molecule_id most be specified')


        experiment = self.client.add_experiment(experiment_machine_id = self.machine_id,
                                                metadata = metadata,
                                                notes = notes,
                                                raw_data_path = raw_data_path,
                                                synthesis_id = synthesis_id,
                                                molecule_id = molecule_id,
                                                parent_experiment_id= parent_experiment_id)                       

        # if synthesis_id:
        #     df = self.client.get('synthesis', filters=[self.client.models.synthesis.synthesis_id == synthesis_id])
        #     molecule_id = df.at[0, 'molecule_id']

        # self.update('experiment', {'molecule_id' : molecule_id}, experiment.uuid)

        return experiment


    def _plot_tile(self, data, title = None, col_num = 5, figsize = None, bar_plot = False, save_filename = None):
        import math
        from matplotlib import gridspec 


        row_num = math.ceil(len(data) / col_num)
        if not figsize:
            figsize = (col_num * 4, row_num * 2.5)

        figs1 =plt.figure(figsize = figsize)
        gs = gridspec.GridSpec(row_num ,col_num)
        #print(data['y_units_id'][0])
        ax = [[] for i in range(len(data))]
        for i, index in enumerate(data.index[::-1]):
            ax[i] = figs1.add_subplot(gs[i //col_num, i % col_num], \
                    xlabel = self.get_unit_name(data['x_units_id'][index]), ylabel = self.get_unit_name(data['y_units_id'][index]))
            if bar_plot:
                ax[i].bar(data['x'][index], data['y'][index])
            else:
                ax[i].plot(data['x'][index], data['y'][index])
            ax[i].legend(['%s' %data['notes'][index]], loc = 'upper right')
            #ax[i].set_title('%s' %data['notes'][index])
        figs1.tight_layout()
        if title:
            figs1.suptitle(title)
            plt.subplots_adjust(top=0.95)
        if save_filename:
            figs1.savefig(save_filename)
            print('%s was saved' %save_filename)
        plt.show()


###device specific functions#####################


    def _parse_metadata(self, metadata):

        data = {
            'filename' : metadata['filename'],
            'Experiment timestamp' : metadata['experiment_timestamp'],
            'excitation_wavelength' : metadata['laser']['wavelength'],
            'excitation_frequency' : metadata['excitation_frequency'],
            'solvent' : metadata['sample']['solvent'],
            'concentration' : metadata['sample']['concentration'],
            'abs_coeff_at_excitation' : metadata['sample']['abs_coeff_at_excitation']
        }

        return data


    def add_TE_results(self, results, synthesis_id = None, molecule_id = None, experiment_id = None, parent_experiment_id = None):
        """
        add transient emission measurement result to the database. Following will be added to the DB.
            - experiment (experiment, if not experiment_id specified)
            - chromatogram (xy_data, with peak_table as metadata)
            - absopriton spectra at the peaks detected in the chromatogram (xy_data)
            - MS spectra at the peaks detected in the chromatogram (xy_data, if result contains the MS data)

        :param result(dict): result dict of the TE measurement (see transient_emission.py)
        :param synthesis_id(uuid or None): id of the synthesis which the experiment to add associated with.
        :param molecule_id(uuid or None): id of the molecule targeted in the experiment.
        :param experiment_id(uuid or None): specify experiment_id in order to add data to the existing experiments.
        :param parants_experiment_id(uuid or None): id of the experiment which the experiment to add associated with.

        returns: experiment_id(uuid) of the added 
        """
        if sum([synthesis_id is not None, molecule_id is not None, experiment_id is not None]) != 1:
            raise Exception('synthesis_id xor molecule_id xor experiment_id should be specified.')
        
        if experiment_id == None:

            metadata = self._parse_metadata(results['metadata'])
            experiment = self.add_experiment(metadata, metadata['filename'], synthesis_id = synthesis_id, \
                                            molecule_id = molecule_id, parent_experiment_id = parent_experiment_id)
            experiment_id = experiment.uuid

            self.update('experiment', {'experiment_timestamp' : metadata['Experiment timestamp']}, experiment.uuid)

        
        #add raw data
        xy_data = self.add_xy_data(
            experiment_id = experiment_id,
            data_name = 'raw data',
            data = [results['raw_data'][0], results['raw_data'][1]],
            x_units_id = self.unit_ids['time'],
            y_units_id = self.unit_ids['signal'],
            metadata = None,
            notes = None)

        #add fitting results
        xy_data = self.add_xy_data(
            experiment_id = experiment_id,
            data_name = 'fitting results',
            data = [results['fitted_data'][0], results['fitted_data'][1]],
            x_units_id = self.unit_ids['time'],
            y_units_id = self.unit_ids['signal'],
            metadata = results['fitting_results'],
            notes = None)
        
        return experiment_id


###some utilities to interact with the DB#########################

    def get_fitting_results(self, experiment_id):
        """
        get fitting result of the experiment

        :param experiment_id (uuid):

        :returns: dict
        """
        df = self.get('xy_data', filters = {self.client.models.xy_data.name == 'fitting results',\
                                                 self.client.models.xy_data.experiment_id == experiment_id})
        if len(df) == 0:
            raise Exception('No fitting results is stored in experiment <%s>' %experiment_id)

        return df['metadata'][0]
 

    def get_raw_data(self, experiment_id):
        """
        get raw data of the experiment

        :param experiment_id (uuid):

        :returns: list
        """
        df = self.get('xy_data', filters = {self.client.models.xy_data.name == 'raw data',\
                                                self.client.models.xy_data.experiment_id == experiment_id})
        if len(df) == 0:
            raise Exception('No raw_data is stored in experiment <%s>' %experiment_id)

        return [df['x'][0], df['y'][0]]


    def get_fitted_data(self, experiment_id):
        """
        get fitted data of the experiment

        :param experiment_id (uuid):

        :returns: list
        """
        df = self.get('xy_data', filters = {self.client.models.xy_data.name == 'fitting results',\
                                                self.client.models.xy_data.experiment_id == experiment_id})
        if len(df) == 0:
            raise Exception('No fitting result is stored in experiment <%s>' %experiment_id)

        return [df['x'][0], df['y'][0]]

 
    def plot_result(self, experiment_id, save_filename = None):
        """
        plot result of the experiment. 

        :param experiment_id (str):
        :param save_filename(str)
        """
        metadata = self.get_experiment(experiment_id = experiment_id)['metadata'][0]
        raw = self.get_raw_data(experiment_id)
        fitted = self.get_fitted_data(experiment_id)
        params = self.get_fitting_results(experiment_id)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, xlim = [10,0.9*1e9/metadata['excitation_frequency']], \
            ylim = [1e0, np.max(raw[1])*2], yscale = 'log', xlabel = 'time/ns', ylabel = 'counts')
        ax.plot(raw[0], raw[1])
        if len(fitted) > 0:
            if 'tau' in params.keys():
                ax.plot(fitted[0], fitted[1])
                ax.text(0.7, 0.85, 'tau = {:.2f} ns'.format(params['tau']), transform=ax.transAxes)
                ax.text(0.7, 0.77, 'R2 = {:.4f}'.format(params['R2']), transform=ax.transAxes)
            if 'tau1' in params.keys():
                ax.plot(fitted[0], fitted[1])
                ax.text(0.7, 0.9, 'tau1 = {:.2f} ns'.format(params['tau1']), transform=ax.transAxes)
                ax.text(0.7, 0.825, 'amp1 = {:.4f}'.format(params['amp1']), transform=ax.transAxes)
                ax.text(0.7, 0.75, 'tau2 = {:.2f} ns'.format(params['tau2']), transform=ax.transAxes)
                ax.text(0.7, 0.675, 'amp2 = {:.4f}'.format(params['amp2']), transform=ax.transAxes)
                ax.text(0.7, 0.6, 'R2 = {:.4f}'.format(params['R2']), transform=ax.transAxes)

        if save_filename:
            plt.savefig(save_filename)
        plt.show()
