from ..Utils import load_json, save_json, timestamp_datetime
from typing import Union, Tuple, List
from pathlib import Path
from .Synthesis import Synthesis
from ..ChemSpeedModules import ISynth
from .Operation import Operation, Communication
from .StepManager import StepManager
import os
import glob
import pathlib


class Batch(object):
    """
    Class to describe a batch of reactions, following a specific synthesis protocol.
    Handles the following steps:
        - Assignment of individual reactions to reaction vessels
        - Preparation of all operations required to complete the synthesis of the respective batch.
        - Iterator to return and execute the next operation as necessary.
    """
    def __init__(self, synthesis_procedure: Path, synthesis_folder: Path, thread="Synthesis"):
        self.synthesis_path: Path = synthesis_folder
        self.synthesis: Synthesis = Synthesis(synthesis_procedure, thread)
        self.thread = thread
        self.batch, self.batch_name, self.batch_file = self._search_batch()
        self.operations = dict()
        self.operations_iterator = None

    def _search_batch(self) -> Tuple[Union[dict, None], Union[str, None], Union[Path, None]]:
        """
        Iterates through all json files in the self.batch_folder directory.
        Loads a json file from that folder, if exists. Returns the loaded dictionary and the file name.
        """
        batch_folder = os.path.join(self.synthesis_path,"Batches_to_Make")
        file_list = list(pathlib.Path(batch_folder).glob('*.json'))
        if file_list:
            return self._load_batch(file_list[0])
        else:
            return None, None, None

    @staticmethod
    def _load_batch(batch_file: Path) -> Tuple[Union[dict, None], Union[str, None], Union[Path, None]]:
        """
        Loads a batch from the corresponding json file.
        Returns the dictionary and the file name.
        """
        batch = load_json(batch_file)
        #batch_name = f"{batch_file.stem}_{timestamp_datetime()}"
        batch_name = f"{batch_file.stem}"
        return batch, batch_name, batch_file

    def _save_batch(self) -> None:
        """
        Saves the current batch including all assigned vials as a json file.
        """
        save_json(self.batch, os.path.join(self.synthesis_path,f"Completed_Batches/{self.batch_name}.json"))

    def _save_status(self) -> None:
        """
        Saves the current status of self.threads as a json file.
        """
        operations_dict = {key: self.operations[key].get_dictionary() for key in self.operations}
        save_json(operations_dict, os.path.join(self.synthesis_path, f"Completed_Batches/{self.batch_name}_Operations.json"))
        
    def clear_batch(self):
        """
        Saves the current (final) status of the batch.
        Deletes the original batch file.
        Clears the current batch.
        """
        self._save_status()
        self.batch_file.unlink()

    def assign_batch(self, reactor: ISynth):
        """
        Assigns all reactions within a batch to a vessel in the reactor object.

        Groups the reactions by occurring variable compounds in order to minimize the amount of transfer steps.
            1. Gets all variable compounds from reaction_ids set, each with a list of reactions that use this compound.
            2. Pick the compound that is used in most reactions.
            3. Assign this group of reactions to a common drawer in the reactor (if possible).
            4. Remove the assigned reactions from reaction_ids
            5. Repeat

        Assigns the positions in the reactor object and in self.batch.
        """
        reaction_ids = set(self.batch.keys())
        reactor.set_batch_vials(len(reaction_ids))

        while len(reaction_ids) > 0:
            groups = self._get_variable_compounds(reaction_ids)
            largest_group_key = max(groups, key=lambda x: len(groups[x]))
            assignment = reactor.assign_group(groups[largest_group_key])
            self._set_location(assignment)
            reaction_ids = reaction_ids - groups[largest_group_key]

        self._save_batch()

    def _get_variable_compounds(self, reactions: set) -> dict:
        """
        Get a dictionary of all variable compounds to complete the synthesis of a given batch.
            - Keys: Compound Identifiers
            - Values: List of all reaction identifiers containing this compound
        """
        compound_types = self.synthesis.get_compound_types()
        all_cmpds = dict()

        for reaction in reactions:
            for compound_type in compound_types:
                cmpd = self.batch[reaction][compound_type]

                if cmpd in all_cmpds:
                    all_cmpds[cmpd].add(reaction)
                else:
                    all_cmpds[cmpd] = {reaction}

        return all_cmpds

    def _set_location(self, locations: dict) -> None:
        """
        Sets the location of a reaction in self.batch according to the specified location in the argument locations.
        """
        for reaction in locations:
            self.batch[reaction]["vial"] = locations[reaction]

    def process_batch(self, reactor: ISynth, step_settings: Path) -> None:
        """
        Parses the batch to eventually get a list of single-step threads to execute sequentially.
            1. Assigns experiments in batch to individual reactors in the reactor object.
            2. Iterates through all steps of the synthesis procedure
                1. Gets a global pre-operation, if applicable
                2. Gets a list of threads for each step
                3. Processes experiments by the StepManager (grouping if appropriate, local pre- and post-threads)
                4. Appends to self.threads
        """

        step_manager = StepManager(step_settings)
        self.assign_batch(reactor)

        self._append_to_operations(Communication(f"Now starting the synthesis of batch {self.batch_name}", self.thread))

        for step in self.synthesis.steps:
            if self.synthesis.perform_global_pre_operation(step):
                pre_operation = step_manager.global_pre_operation(**self.synthesis.get_step_details(step))
                self._append_to_operations(*pre_operation)

            operations = self.synthesis.get_operations(step, self.batch)
            operations = step_manager.process_operations(operations, reactor, self.synthesis.step_is_local(step))
            self._append_to_operations(*operations)

        self._append_to_operations(Communication(f"Synthesis of batch {self.batch_name} has been completed", self.thread))

        self.operations_iterator = iter(self.operations)

    def _append_to_operations(self, *args) -> None:
        """
        Appends all threads sequentially to self.threads.
        """
        index = len(self.operations.keys())
        for operation in args:
            index += 1
            self.operations[index] = operation

    def run_next_operation(self) -> Tuple[str, list, str, dict]:
        """
        Get the next operation (task and parameters) from the self.threads iterator.
        Runs the execute function on the operation to set a timestamp for execution.
        Saves the updated self.threads dictionaroy.
        """
        next_index = next(self.operations_iterator)
        task, target, source, parameters = self.operations[next_index].execute()
        self._save_status()
        return task, target, source, parameters

    def get_operations_list(self) -> List[Operation]:
        """
        Returns a list of all threads (i.e. of all values in self.threads).
        """
        return list(self.operations.values())
