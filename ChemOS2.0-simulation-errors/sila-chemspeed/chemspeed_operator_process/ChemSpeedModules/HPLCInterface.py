from .ChemSpeedModules import ChemSpeedModule
from pathlib import Path
from ..Utils.FileHandling import load_pkl, save_pkl
import time
from typing import Union, Iterable, List
from sila2.client import SilaClient


class HPLCInterface(ChemSpeedModule):
    """
    Interface with the ThermoFisher HPLC-MS and the respective control code by Kazu.
        - instrument status to be read from pkl file
        - job submission via writing of pkl files
        - job results via reading of pkl files
    """
    def __init__(self, name: str, role: str, output_path: Path, positions: Iterable, output_folder: str, defaults_path: Path, **kwargs):
        """
        Instantiates the HPLC interface by calling the __init__ method of the ChemSpeedModule MetaClass
        Additionally sets the path attributes.
            - self.status_file (pkl file containing the current HPLC status)
            - self.submission_path (path to the folder for job submissions)
            - self.result_path (path to the folder for calculation results)
            - self.injection_port
            - self.wash_port
        """
        super().__init__(name, role, output_path, positions, output_folder, defaults_path)
        self.status_file = Path(kwargs["status_file"])
        self.submission_path = Path(kwargs["submission_path"])
        self.result_path = Path(kwargs["result_path"])
        self.injection_port = self.positions
        self.wash_port = kwargs["wash_port"]

        print(f"RESULT PATH IS", self.result_path)

        print(f"HPLC STATUS IS", self._get_hplc_status())

    def wait_for_status(self, status: str) -> None:
        """
        Waits until the HPLC-MS instrument is in a given status.
        """
        hplc_status = "socket connection offline"

        return  hplc_status

        while hplc_status != status:
            hplc_status = self._get_hplc_status()
            print(f"HPLC STATUS IS: {hplc_status}")
        return

    def _get_hplc_status(self) -> str:
        """
        Gets the current status of the HPLC-MS instrument, as specified in the status file.
        """

        return "socket connection offline"
        try:
            client = SilaClient("127.0.0.1", 65004, insecure=True)
            status = client.HPLCMSsimulator.Status().Termination
            return status
        except (FileNotFoundError, KeyError, PermissionError):
            return "fail"

    def get_hplc_results(self, job_name: str) -> Union[dict, None]:
        """
        Reads the output of a HPLC-MS job. Returns the corresponding dictionary.
        """
        # for _ in range(30):
        #     try:
        #         hplc_results = load_pkl(self.result_path / f"{job_name}.pkl")
        #         print(f"DEBUG: HPLC Results were read in: {hplc_results}")
        #         return hplc_results
        #     except (FileNotFoundError, PermissionError):
        #         time.sleep(1)

        return None

    def submit_hplc_job(self, experiment_identifier: str, injection_type: str, targets: List[dict], **kwargs) -> None:
        """
        Writes a pickle file to submit a job for HPLC-MS analysis.
        """
        job_settings = {
            "name": f"{experiment_identifier}_{injection_type}",
            "itype": injection_type,
            "targets": targets
        }
        print(f"DEBUG: HPLC Job Submitted: {job_settings}")
        save_pkl(job_settings, self.submission_path / f"{experiment_identifier}_{injection_type}.pkl")

        # def printupdate(Submit_Job_Response) -> None:
        #     print(Submit_Job_Response.Status)
        # client = SilaClient("127.0.0.1",65004, insecure=True)
        # instance = client.HPLCMSsimulator.SubmitJob("procedure.txt")
        # sub = instance.subscribe_to_intermediate_responses()
        # sub.add_callback(printupdate)
        # while instance.done != True:
        #     time.sleep(1)
        # try:
        #     if instance.done == True:
        #         print("job completed")    
        # except:
        #     print("error logging final info")
