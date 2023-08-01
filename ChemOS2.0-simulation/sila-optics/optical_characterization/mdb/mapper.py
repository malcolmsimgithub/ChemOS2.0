import logging

from . import utils


def schema_decorator(func):
    # sqlalchemy.exc.IntegrityError: (psycopg2.errors.UniqueViolation)
    pass


class SchemaMapper:
    """
    SchemaMapper map the schema of the database. It contains some utility function to
    add data to the database.

    It is meant to be used through MDBClient.
    """
    def __init__(self, dao, logger=None):
        self.dao = dao
        self.logger = logger or logging.getLogger(__name__)

    def add_molecule_type(self, name):
        data = {'name': name}
        event = self.dao.add('molecule_type', data)
        f_event = self.dao.commit_or_fetch_event('molecule_type', data)
        if f_event:
            self.logger.warning(f"molecule_type {name} already exists (event: {f_event.uuid})")
            return f_event
        return event

    def add_molecule(self, smiles, molecule_type_id, reactant_id=[],
                     pubchem_autofill=True):
        data = {'smiles': smiles,
                'molecule_type_id': molecule_type_id}

        if pubchem_autofill is True:
            pubchem_data = utils.pubchem_lookup(smiles)
            if pubchem_data is not None:
                data.update(pubchem_data)

        event = self.dao.add('molecule', data)
        f_event = self.dao.commit_or_fetch_event('molecule', data)
        if f_event:
            self.logger.warning(f"Molecule {smiles} already exists (event: {f_event.uuid})")
            return f_event

        if not isinstance(reactant_id, list):
            reactant_id = [reactant_id]
        for order, id in enumerate(reactant_id):
            data = {'product_molecule_id': event.uuid,
                    'reactant_molecule_id': id,
                    'order': order}
            self.dao.add('molecule_molecule', data)
        self.dao.session.commit()
        return event

    def add_conformer(self, molecule_id, x, y, z, atomic_numbers, metadata):
        data = {'x': x, 'y': y, 'z': z, 'atomic_numbers': atomic_numbers, 'metadata': metadata}
        event = self.dao.add('conformer', data)
        self.dao.session.commit()
        if not isinstance(molecule_id, list):
            molecule_id = [molecule_id]
        for m_id in molecule_id:
            data = {'molecule_id': m_id,
                    'conformer_id': event.uuid}
            self.dao.add('conformer_molecule', data)
        self.dao.session.commit()
        return event

    def add_calculation_type(self, name):
        data = {'name': name}
        event = self.dao.add('calculation_type', data)
        f_event = self.dao.commit_or_fetch_event('calculation_type', data)
        if f_event:
            self.logger.warning(f"calculation_type {name} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_calculation(self, input, output, command_line, calculation_type_id,
                        software_id, conformer_id, metadata, output_conformer_id=None):
        data = {'input': input,
                'output': output,
                'command_line': command_line,
                'calculation_type_id': calculation_type_id,
                'software_id': software_id,
                'conformer_id': conformer_id,
                'metadata': metadata}
        if output_conformer_id is not None:
            data['output_conformer_id'] = output_conformer_id
        event = self.dao.add('calculation', data)
        self.dao.session.commit()
        return event

    def add_software(self, name, version):
        data = {'name': name, 'version': version}
        event = self.dao.add('software', data)
        f_event = self.dao.commit_or_fetch_event('software', data)
        if f_event:
            self.logger.warning(f"sofware {name}-{version} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_lab(self, name, short_name):
        data = {'name': name, 'short_name': short_name}
        event = self.dao.add('lab', data)
        f_event = self.dao.commit_or_fetch_event('lab', data)
        if f_event:
            self.logger.warning(f"lab {name} - {short_name} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_synthesis_machine(self, name, make, model, metadata, lab_id):
        data = {'name': name,
                'metadata': metadata,
                'lab_id': lab_id,
                'make': make,
                'model': model}
        event = self.dao.add('synthesis_machine', data)
        f_event = self.dao.commit_or_fetch_event('synthesis_machine', data)
        if f_event:
            self.logger.warning(f"synthesis_machine {name}-{make}-{model} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_synthesis(self, synthesis_machine_id, targeted_molecule_id, xdl, notes):
        data = {'molecule_id': targeted_molecule_id,
                'synthesis_machine_id': synthesis_machine_id,
                'xdl': xdl,
                'notes': notes}
        event = self.dao.add('synthesis', data)
        self.dao.session.commit()
        return event

    def add_synthesis_molecule(self, synthesis_id, molecule_id, yield_):
        data = {'synthesis_id': synthesis_id,
                'molecule_id': molecule_id,
                'yield': yield_}
        event = self.dao.add('synth_molecule', data)
        self.dao.session.commit()
        return event

    def add_experiment_machine(self, name, make, model, metadata, experiment_type_id, lab_id):
        data = {'name': name,
                'lab_id': lab_id,
                'experiment_type_id': experiment_type_id,
                'make': make,
                'model': model,
                'metadata': metadata}
        event = self.dao.add('experiment_machine', data)
        f_event = self.dao.commit_or_fetch_event('experiment_machine', data)
        if f_event:
            self.logger.warning(f"experiment_machine {name}-{make}-{model} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_experiment_type(self, name):
        data = {'name': name}
        event = self.dao.add('experiment_type', data)
        f_event = self.dao.commit_or_fetch_event('experiment_type', data)
        if f_event:
            self.logger.warning(f"experiment_type {name} already exist (event: {f_event.uuid})")
            return f_event
        return event

    def add_experiment(self, experiment_machine_id, metadata,
                       notes, raw_data_path, synthesis_id=None,
                       molecule_id=None, parent_experiment_id=None):
        data = {'experiment_machine_id': experiment_machine_id,
                'metadata': metadata,
                'raw_data_path': raw_data_path,
                'notes': notes}
        if (synthesis_id is None and molecule_id is None) or \
           (synthesis_id is not None and molecule_id is not None):
            raise ValueError('Either synthesis_id xor molecule_id most be specified')
        if synthesis_id is not None:
            data['synthesis_id'] = synthesis_id
        if molecule_id is not None:
            data['molecule_id'] = molecule_id
        if parent_experiment_id is not None:
            data['parent_experiment_id'] = parent_experiment_id
        event = self.dao.add('experiment', data)
        self.dao.session.commit()
        return event

    def add_xy_data_experiment(self, experiment_id, name, x, y, x_units_id, y_units_id):
        data = {'experiment_id': experiment_id,
                'name': name,
                'x': x,
                'y': y,
                'x_units_id': x_units_id,
                'y_units_id': y_units_id
                }
        event = self.dao.add('xy_data', data)
        self.dao.session.commit()
        return event

    def add_xy_data_calculation(self, calculation_id, name, x, y, x_units_id, y_units_id):
        data = {'calculation_id': calculation_id,
                'name': name,
                'x': x,
                'y': y,
                'x_units_id': x_units_id,
                'y_units_id': y_units_id}
        event = self.dao.add('xy_data', data)
        self.dao.session.commit()
        return event

    def add_xyz_data_experiment(self, experiment_id, name, x, y, z,
                               x_units_id, y_units_id, z_units_id):
        data = {'experiment_id': experiment_id,
                'name': name,
                'x': x,
                'y': y,
                'z': z,
                'x_units_id': x_units_id,
                'y_units_id': y_units_id,
                'z_units_id': z_units_id}
        event = self.dao.add('xyz_data', data)
        self.dao.session.commit()
        return event

    def add_xyz_data_calculation(self, calculation_id, name, x, y, z,
                                x_units_id, y_units_id, z_units_id):
        data = {'calculation_id': calculation_id,
                'name': name,
                'x': x,
                'y': y,
                'z': z,
                'x_units_id': x_units_id,
                'y_units_id': y_units_id,
                'z_units_id': z_units_id}
        event = self.dao.add('xyz_data', data)
        self.dao.session.commit()
        return event

    def add_data_unit(self, name):
        data = {'name': name}
        event = self.dao.add('data_unit', data)
        f_event = self.dao.commit_or_fetch_event('data_unit', data)
        if f_event:
            self.logger.warning(f"data_unit {name} already exist (event: {f_event.uuid})")
            return f_event
        return event
