from ..Workflows import Operation
from ..Utils import save_csv, load_json, save_json, timestamp_datetime
from ..ChemSpeedModules import *
import pandas as pd
from pathlib import Path
import os
import json


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
                 hardware: dict
                 ):
        """
        Instantiates the HPLCCharacterizationScheduler.

        Parameters:
            default_settings: Path to the folder where default_settings objects are stored.
            characterization_path: Path to the folder where external files for characterization are saved.
            hardware: Dictionary of instrument hardware

        Sets the following attributes:
            self._settings (dict): Dictionary of all HPLC characterization-related settings.
            self.files (dict): Dictionary of all regularly written output files.
            self.hplc_interface
            self.storage_rack
            self.filtration_rack: Hardware components required for operating the HPLCCharacterizationScheduler.
            self._injection_counter (int): Counter for regular blank injections.
            self.characterizations (dict): Dictionary of all previous and current characterization runs.
            self.operations (list): List of operations to be performed by the executor.
            self._running (str): Identifier of the currently _running characterization.
        """

        self.characterization_path = characterization_path
        self._settings: dict = load_json(os.path.join(default_settings,"HPLCCharacterizationScheduler_Settings.json"))

        self.files: dict = {
            "characterizations": os.path.join(characterization_path,"Characterizations.json"),
            "jobs": os.path.join(characterization_path,"Positions_for_Characterization.csv"),
            "queue": os.path.join(characterization_path,"Characterization_Queue.csv")
        }
        
        self._output_path: Path = os.path.join(characterization_path,"Completed")

        if not os.path.exists(self._output_path):
            os.makedirs(self._output_path)

        

        self.hplc_interface, self.storage_rack, self.filtration_rack = self._get_hardware_components(hardware)

        self._injection_counter: int = 0

        self.characterizations = dict()
        self._load_characterizations()
        self._generate_job_file()

        self.queue = list()
        self._update_queue()

        self.operations = list()
        self._running = None

        self._update_all_files()

    def _get_hardware_components(self, hardware: dict) -> tuple:
        """
        Checks the entire hardware dictionary and identifies the necessary hardware components by the "role" attribute.
            - hplc_interface: "hplc_interface",
            - storage_rack: "sample_storage",
            - filtration_rack: "filtration_rack"
        """
        hplc_interface, storage_rack, filtration_rack = None, None, None
        for component in hardware.values():
            if component is not None:
                if component.role == "analysis":
                    hplc_interface = component
                elif component.role == "sample":
                    storage_rack = component
                elif component.role == "filtration":
                    filtration_rack = component

        if hplc_interface and storage_rack and filtration_rack:
            return hplc_interface, storage_rack, filtration_rack
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
            self._load_next_characterization()
            return self.execute_next_operation()

    # Methods for internal and external file handling (queue, characterizations, job file)

    def _save_queue_file(self):
        """
        Saves the current queue file.
        """
        queue_for_writing = {str(i): (reaction, self.characterizations[reaction]["structure"]) for i, reaction in enumerate(self.queue)}
        save_csv(queue_for_writing, self.files["queue"], ("Number", "Identifier", "Structure"))

    def _load_characterizations(self) -> None:
        """
        Loads the previously saved characterization dictionary from the corresponding json file.
        """
        if os.path.exists(self.files["characterizations"]):
            self.characterizations = load_json(self.files["characterizations"])

    def _save_characterization(self) -> None:
        """
        Saves the current self.characterizations dictionary to the corresponding json file.
        """
        #save_json(self.characterizations, self.files["characterizations"])

    def _update_all_files(self):
        """
        Loads / saves / updates all files
        """
        self._update_submitted_jobs()
        self._update_queue()
        self._save_characterization()
        self._save_queue_file()

    def _generate_job_file(self) -> None:
        """
        Method for starting up the CharacterizationScheduler and setting up the self.files["jobs"] file.
        Iterates over self.characterizations and re-assigns all storage vials to the respective compounds.
        Creates a pd.DataFrame and saves it as a csv file.
        """
        internal_df = pd.DataFrame([], columns=["Identifier", "Structure", *self._settings["operations"].keys()])

        for vial in self.storage_rack.get_positions():

            # Checking if a non-completed characterization exists for this vial
            identifier = None
            for characterization in self.characterizations:
                if self.characterizations[characterization]["vial"] == vial and self.characterizations[characterization]["status"] != "completed":
                    self.storage_rack.update_storage(vial, characterization)
                    identifier = characterization
                    break

            # Writing the characterization details or an empty line for each vial
            if identifier:
                internal_df.loc[vial] = [
                    identifier,
                    self.characterizations[identifier]["structure"],
                    *[self.characterizations[identifier][key] for key in self._settings["operations"]]
                ]
            else:
                internal_df.loc[vial] = ["" for _ in internal_df.columns]

        internal_df.to_csv(self.files["jobs"], index_label="Position")

    def _generate_output_file(self, job: str) -> None:
        """
        Generates an output file that states the completion of the job.

        Args:
            job: Name of the job
        """
        output: dict = {
            "name": job,
            "identifier": self.characterizations[job]["identifier"],
            "smiles": self.characterizations[job]["structure"],
            "synthesis_result": self.characterizations[job]["detected"]
        }

        with open(os.path.join(self.characterization_path,f"Completed/{job}.json"), 'w') as jsonfile:
            json.dump(output, jsonfile)
        


    # Methods to deal with characterization scheduling (loading, queuing, ...)

    def _update_submitted_jobs(self):
        """
        Manages the submitted jobs csv file.
            - adds new characterization jobs to self.characterizations and self.queue
        """
        loaded_file = pd.read_csv(self.files["jobs"], index_col=0)

        for vial in loaded_file.index.values:
            if pd.isna(loaded_file.at[vial, "Identifier"]):
                continue

            exp_name = loaded_file.at[vial, "Identifier"]

            if exp_name in self.characterizations:
                continue

            elif exp_name != "None":
                self._add_characterization(
                    exp_name=exp_name,
                    identifier=loaded_file.at[vial, "Identifier"],
                    vial=vial,
                    structure=loaded_file.at[vial, "Structure"],
                    # TODO: Implement loading of multiple target compounds? -> Photocatalysis?
                    **{key: bool(loaded_file.at[vial, key]) for key in self._settings["operations"]}
                )

    def _add_characterization(self, exp_name: str, **kwargs) -> None:
        """
        Adds new characterization to be run.
        If the experiment identifier is already present in self.characterizations, returns None immediately.

        Updates self.characterizations, self.storage_rack, self.queue

        Parameters:
            exp_name (str): identifier of the characterization run (id_timestamp)
            **kwargs: keyword arguments to be added to the dictionary in self.characterizations[identifier]
                        - vial (str)
                        - structure (str)
                        - $OPERATION (bool)
        """
        if exp_name in self.characterizations:
            return

        self.characterizations[exp_name] = {
            "status": "queuing",
            "detected": False,
            "completed": {}
        }
        self.characterizations[exp_name].update(kwargs)

        self.storage_rack.update_storage(kwargs["vial"], exp_name)
        self._update_queue()

    def _update_queue(self) -> None:
        """
        Iterates over all characterizations in self.characterizations and adds to self.queue if
            - not queued already
            - status is "queuing"
        """
        for characterization in self.characterizations:
            if characterization not in self.queue and self.characterizations[characterization]["status"] == "queuing":
                self.queue.append(characterization)

    def _load_next_characterization(self) -> None:
        """
        Loads the next characterization into self._running.
        Raises a StopIteration if no further characterizations are found in the queue.
        """
        if self.queue:
            self._running = self.queue.pop(0)
            self.characterizations[self._running]["jobs"] = [job for job in self._settings["operations"] if self.characterizations[self._running][job]]

        else:
            self._update_submitted_jobs()
            if self.queue:
                self._load_next_characterization()
            else:
                raise StopIteration("No characterizations found!")

    def _get_new_operations(self) -> None:
        """
        Get the next operations, as specified in self._running["jobs"], and append them to self.operations.
        Raises a ValueError if self._running["jobs"] is empty.
        """
        if self._running is None:
            raise ValueError("self._running is currently empty!")
        if not self.characterizations[self._running]["jobs"]:
            self._complete_characterization()
            raise ValueError("No operations left in self._running!")

        job: str = self.characterizations[self._running]["jobs"][0]

        # TODO: Think of a good way to make these hard-coded operations a bit more flexible

        if job == "filter_collect":
            try:
                self._filter_collect_operations()
            except PositionNotAvailableError:
                return

        if job == "characterization_1st" or "characterization_2nd":
            self.hplc_interface.wait_for_status("idle")
            self._hplc_inject_operation(job)

        # if job == "characterization_2nd":
        #     self.hplc_interface.wait_for_status("idle")
        #     if self._check_for_product_detection(self._running):
        #         self.characterizations[self._running]["detected"] = True
        #         self._hplc_inject_operation(job)
        #     else:
        #         communication = Operation(
        #             experiment_identifier=[self._running],
        #             task="communicate",
        #             parameters={
        #                 "thread": "Characterization",
        #                 "message": f"Injection 'characterization_2nd' skipped for {self._running}"
        #             }
        #         )
        #         self._append_to_operations(communication)

        self.characterizations[self._running]["completed"][job] = timestamp_datetime()
        del self.characterizations[self._running]["jobs"][0]

    def _complete_characterization(self) -> None:
        """
        Complete a characterization run by updating self.characterizations and setting self._running back to None.
        Saves the current status of self.characterizations.
        """
        self.characterizations[self._running]["status"] = "completed"

        # Clear the entry in the submission file
        submission_file = pd.read_csv(self.files["jobs"], index_col=0)
        for column in submission_file.columns.values:
            submission_file.at[self.characterizations[self._running]["vial"], column] = ""
        submission_file.to_csv(self.files["jobs"])

        # Move entry to new key with timestamp of completion
        # new_id: str = f"{self._running}_{timestamp_datetime()}"
        new_id: str = f"{self._running}"
        self.characterizations[new_id] = self.characterizations.pop(self._running)
        #self._generate_output_file(new_id)

        # Communicate completion
        communication = Operation(
            experiment_identifier=[self._running],
            task="communicate",
            parameters={
                "thread": "Characterization",
                "message": f"Characterization of {self._running} completed."
            }
        )
        self._append_to_operations(communication)
        self._running = None

        self._update_all_files()

    def _append_to_operations(self, *args, insert_at_beginning=False) -> None:
        """
        Appends an operation to self.operations.
        """
        if insert_at_beginning:
            self.operations = [*args] + self.operations
        else:
            self.operations = self.operations + [*args]

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
        if "filtration" not in self.characterizations[self._running]:
            filtration_position = self.filtration_rack.get_filtration_position(self._running)
            self.characterizations[self._running]["filtration"] = filtration_position

    def _filter_collect_operations(self):
        """
        Generates the Operation object required for performing a filter_collect step.

        If no position is available on the FiltrationRack (PositionNotAvailableError),
        a cleaning confirmation request is required.
        """
        try:
            self._set_filtration_position()
            filter_collect = Operation(
                experiment_identifier=[self._running],
                source=self.characterizations[self._running]["vial"],
                target=self.characterizations[self._running]["filtration"],
                **self._settings["operations"]["filter_collect"]
            )
            self._append_to_operations(filter_collect)

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

    def _hplc_inject_operation(self, injection_type: str) -> None:
        """
        Generates the Operation objects required for performing an HPLC injection.
            - submit_job
            - inject_to_hplc (either from self._running["vial"] or from self._running["filtration"])
        """
        communication = Operation(
            experiment_identifier=[self._running],
            task="communicate",
            parameters={
                "thread": "Characterization",
                "message": f"Injection {injection_type} started for {self._running}"
            }
        )

        # submit_file = Operation(
        #     experiment_identifier=[self._running],
        #     task="submit_hplc_job",
        #     parameters={
        #         "experiment_identifier": self._running,
        #         "injection_type": injection_type,
        #         "targets": [
        #             {
        #                 "name": self.characterizations[self._running]["identifier"],
        #                 "type": "product",
        #                 "smiles": self.characterizations[self._running]["structure"]
        #             }
        #         ]
        #     }
        # )

        if "filtration" in self.characterizations[self._running]:
            injection = Operation(
                experiment_identifier=[self._running],
                task="inject_to_hplc",
                parameters={
                    "filtration": self.characterizations[self._running]["filtration"]
                }
            )
        else:
            injection = Operation(
                experiment_identifier=[self._running],
                task="inject_to_hplc",
                source=self.characterizations[self._running]["vial"],
            )

        self._append_to_operations(communication, injection)

    def _check_for_product_detection(self, job_name: str) -> bool:
        """
        Evaluates a previous HPLC run (from "characterization_1st") to check if a target_zone compound was found.
        """
        hplc_results = self.hplc_interface.get_hplc_results(f"{job_name}_characterization_1st")

        if not hplc_results:
            return False

        if "concentration" in hplc_results:
            if hplc_results["concentration"]:
                return True
        else:
            return False
