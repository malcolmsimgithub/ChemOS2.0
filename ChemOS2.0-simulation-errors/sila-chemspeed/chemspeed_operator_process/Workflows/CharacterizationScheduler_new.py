from ..Workflows import Operation
from ..Utils import save_csv, load_json, save_json, timestamp_datetime
from ..ChemSpeedModules import *
import pandas as pd
from pathlib import Path
import os
import datetime
import json
import shutil


class HPLCCharacterizationScheduler:
    """
    HPLCCharacterizationScheduler
    Object for flexibly handling (internally stored and external) sample preparation and injection to HPLC-MS.
        - loading of samples from an externally editable csv file ("naive" user interface)
        - queuing and job scheduling for sample preparation and HPLC injection(s)
        - management of Operation objects for higher-level controller (via execute_next_operation) method
    """
    def __init__(self,
                 default_settings: Path,
                 characterization_path: Path,
                 output_path: Path,
                 hardware: dict
                 ):
        """
        Instantiates the HPLCCharacterizationScheduler.

        Parameters:
            default_settings: Path to the folder where default_settings objects are stored.
            characterization_path: Path to the folder where external files for characterization are saved.
            hardware: Dictionary of instrument hardware

        Sets the following attributes:
            self._settings (dict): Dictionary of all HPLC-injection and filtration related settings.
            self.files (dict): Dictionary of all regularly written output files.
            self.storage_rack
            self.filtration_rack: Hardware components required for operating the HPLCCharacterizationScheduler.
            self._injection_counter (int): Counter for regular blank injections.
            self.operations (list): List of operations to be performed by the executor.
            self._running (str): Identifier of the currently _running job.
            self.output_path(Path): output path of jobs
        """
        self.output_path = output_path
        self.characterization_path = characterization_path
        self._settings: dict = load_json(default_settings/"HPLCCharacterizationScheduler_Settings.json")

        self.files: dict = {
            "characterizations": characterization_path/"Characterizations.json",
            "jobs": characterization_path/"Positions_for_Characterization.csv",
            "queue": characterization_path/"Characterization_Queue.csv"
        }
        


        self.storage_rack, self.filtration_rack = self._get_hardware_components(hardware)
        self._injection_counter: int = 0
        self.operations = list()
        self._running = None

        # self._update_all_files()

    def _get_hardware_components(self, hardware: dict) -> tuple:
        #removed hplc interface
        """
        Checks the entire hardware dictionary and identifies the necessary hardware components by the "role" attribute.
            - storage_rack: "sample_storage",
            - filtration_rack: "filtration_rack"
        """
        storage_rack, filtration_rack = None, None
        for component in hardware.values():
            if component is not None:
                if component.role == "sample":
                    storage_rack = component
                elif component.role == "filtration":
                    filtration_rack = component

        if storage_rack and filtration_rack:
            return storage_rack, filtration_rack
        else:
            raise UnknownHardwareException()

    # Public method for executing all operations

    def execute_next_operation(self) -> tuple:
        """
        Executes and returns the next scheduled operation to be performed for Characterization.

        Pops self.operations[0] and runs the execute function on that operation.

        If no operations are queued, triggers the generation of the next set of operations from self._running.
        Otherwise, triggers the generation of new operations
            - a) next characterization operation from self._running
            - b) next sample to characterize from self.queue
        """
        self._check_injection_counter()

        if self.operations:
            operation = self.operations.pop(0)
            self._update_injection_counter(operation)
            return operation.execute()

        try:
            self._get_new_operations()
            return self.execute_next_operation()
        except ValueError:
            return self.execute_next_operation()

    # Methods for internal and external file handling (queue, characterizations, job file)
   

    def _generate_output_file(self, job: str) -> None:
        """
        Generates an output file that states the completion of the injection.

        Args:
            job: Name of the job
        """
        output: dict = {
            "name": job,
            "identifier": self.characterizations[job]["identifier"],
        }

        with open(os.path.join(self.characterization_path,f"Completed/{job}.json"), 'w') as jsonfile:
            json.dump(output, jsonfile)
        
    # Methods to deal with characterization scheduling (loading, queuing, ...)

    def _get_new_operations(self) -> None:
        """
            searches the input folder of Characterization/Characterizations_to_do and adds new operations 
            to the queue for a particular job. sorts by job submission date/time
        """

        job_files = [pos_json for pos_json in os.listdir(self.characterization_path) if pos_json.endswith('.json')]
        job_files.sort(key=lambda x: datetime.datetime.strptime(x[-19:-5], '%y-%m-%d_%H-%M'))
        
        

        if len(job_files) == 0:
            raise StopIteration("self._running is currently empty!")
        else:
            newjobfile = job_files[0]
        
        with open(self.characterization_path/newjobfile, "r") as f:
            job = json.load(f)
            job['path'] = self.characterization_path/newjobfile

        
        self._running = job
        
        shutil.move(self.characterization_path/newjobfile, self.output_path/newjobfile)


        # TODO: Think of a good way to make these hard-coded operations a bit more flexible

        if job['Operation'] == "filter_collect":
            try:
                self._filter_collect_operations()
            except PositionNotAvailableError:
                return

        if job['Operation'] == "inject_to_hplc":
            self._hplc_inject_operation()
        




    def _append_to_operations(self, *args, insert_at_beginning=False) -> None:
        """
        Appends an operation to self.operations.
        """
        if insert_at_beginning:
            self.operations = [*args] + self.operations
        else:
            self.operations = self.operations + [*args]
        print("self.operations:")
        print(len(self.operations))

    # Methods to deal with the injection counter & regular blank injections

    def _update_injection_counter(self, operation: Operation) -> None:
        """
        Updates the injection counter to allow for regular blank / cleaning runs by checking if the
        operation is an injection.

        Parameters:
            operation: Operation to be performed
        """
        if operation.task == "inject_to_hplc":
            self._injection_counter += 1

    def _check_injection_counter(self) -> None:
        """
        Checks the injection counter variable.
        If the max. number of sequential injections (from self._settings) is reached, returns a blank injection and
        re-sets the injection counter to 0.
        Returns None otherwise.
        """
        if self._injection_counter >= self._settings["regular_blank_runs"]:
            self._injection_counter = 0
            self._blank_run()

    def _blank_run(self) -> None:
        """
        Creates the operations to submit a blank run for washing the HPLC-MS.
        Inserts them to the first position of the characterization queue.
        """
        timestamp = timestamp_datetime()
        submit_file = Operation(
            experiment_identifier=[f"blank_{timestamp}"],
            task="submit_hplc_job",
            parameters={
                "experiment_identifier": f"blank_{timestamp}",
                "injection_type": "characterization_1st",
                "targets": [{"name": "MeCN", "type": "product", "smiles": "CC#N"}]}
        )
        injection = Operation(
            experiment_identifier=[f"blank_{timestamp}"],
            task="inject_to_hplc",
            source=self._settings["blank"],
            parameters={"filtration": False},
        )
        self._append_to_operations(submit_file, injection, insert_at_beginning=True)

    # Methods that are specific to certain classes of HPLC/characterization-related operations

    def _set_filtration_position(self) -> None:
        """
        Assigns a filtration position to the entry for self._running.
        """
        if "filtration" not in self._running:
            filtration_position = self.filtration_rack.get_filtration_position(self._running["Identifier"])
            self._running["filtration"] = filtration_position

    def _filter_collect_operations(self):
        """
        Generates the Operation object required for performing a filter_collect step.

        If no position is available on the FiltrationRack (PositionNotAvailableError),
        a cleaning confirmation request is required.
        """
        try:
            self._set_filtration_position()
            communication = Operation(
                experiment_identifier=[self._running['Identifier']],
                task="communicate",
                parameters={
                    "thread": "Characterization",
                    "message": f"Filtering compound for {self._running['Identifier']} in position {self._running['filtration'].replace(':', '')}"
                }
            )
            filter_collect = Operation(
                experiment_identifier=[self._running],
                source=self._running["Position"],
                target=self._running["filtration"],
                **self._settings["operations"]["filter_collect"]
            )
            communication_termination = Operation(
                experiment_identifier=[self._running['Identifier']],
                task="communicate",
                parameters={
                    "thread": "Characterization",
                    "message": f"Filtering done for {self._running['Identifier']}"
                }
            )
            self._append_to_operations(communication, filter_collect, communication_termination)

        except PositionNotAvailableError:
            request_operation = Operation(
                experiment_identifier=[None],
                task="request_confirmation",
                parameters={
                    "message": "The filtration rack has no positions available. Please re-set the filtration positions.",
                    "thread": "Characterization"
                }
            )
            self._append_to_operations(request_operation, insert_at_beginning=True)

            raise PositionNotAvailableError("No filtration position is available.")

    def _hplc_inject_operation(self) -> None:
        """
        Generates the Operation objects required for performing an HPLC injection.
            - inject_to_hplc (either from self._running["vial"] or from self._running["filtration"])
        """

        #M.S: removed hplc submit job function 
        communication = Operation(
            experiment_identifier=None,
            task="communicate",
            parameters={
                "thread": "Characterization",
                "message": f"Injection started for job {self._running['Identifier']}"
            }
        )
        injection = Operation(
            experiment_identifier=[self._running],
            task="inject_to_hplc",
            source=self._running["Position"],
        )

        communication_termination = Operation(
                experiment_identifier=[self._running['Identifier']],
                task="communicate",
                parameters={
                    "thread": "Characterization",
                    "message": f"Injection done for {self._running['Identifier']}"
                }
            )

        self._append_to_operations(communication, injection, communication_termination)
