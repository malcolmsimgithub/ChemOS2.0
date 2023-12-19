import time
from pathlib import Path
from typing import Optional
import socket
import shutil
import os
import json
import random
import sys

INPUT = os.path.abspath("communication/jobs_submitted/")
OUTPUT = os.path.abspath("communication/jobs_returned/")
STATUS = os.path.abspath("communication/HPLCMS_status/HPLCMS_status.txt")  
STORAGE = os.path.abspath("communication/storage")  


HEADER_LENGTH = 10
IP = "127.0.0.1"
PORT = 65011

class HPLCSimulator(object):

    def __init__(self, communication_dir: Path):

        # Define some of the folders required for communication with the Chemspeed and create the sub-folders,
        # if required. Not nicely written, but it should do the job.
        self._status_dir: Path = communication_dir / "HPLCMS_status"
        self._status_dir.mkdir(parents=True, exist_ok=True)
        self._status_file: Path = self._status_dir / "HPLCMS_status.txt"
        self._input_dir: Path = communication_dir / "jobs_submitted"
        self._input_dir.mkdir(parents=True, exist_ok=True)
        self._output_dir: Path = communication_dir / "jobs_returned"
        self._output_dir.mkdir(parents=True, exist_ok=True)

        self._instrument_status: str = "offline"
        self._current_job: Optional[dict] = None

        self.client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Connect to a given ip and port
        self.client_socket.connect((IP, PORT))

        # Set connection to non-blocking state, so .recv() call won;t block, just return some exception we'll handle
        self.client_socket.setblocking(False)

        id  = 'HPLCMS'.encode('utf-8')
        username_header = f"{len(id):<{HEADER_LENGTH}}".encode('utf-8')
        self.client_socket.send(username_header + id)

        

    @property
    def _status(self) -> str:
        return self._instrument_status

    @_status.setter
    def _status(self, target_status) -> None:
        self._instrument_status = target_status
        self._write_status()

    def _write_status(self) -> None:

        status_string = f"{self._instrument_status}"

        self.silacommunicate(f"hplc is: {status_string}")

        with open(self._status_file, "w") as f:
            f.write(status_string)

    def operate(self):

        self._status = "idle"

        while True:


            self._wait_for_input()

            if self._current_job['name'].startswith("blank_run"):
                self.silacommunicate("performing blank run")
                self.blank_run()
                
            else:
                self.silacommunicate("beginning HPLC characterization")
                self._prepare_run()
                self.silacommunicate("waiting for injection from chemspeed...")      
                self._hplc_run()
                self._data_analysis()
            

    def _wait_for_input(self) -> None:
        """
            Constantly checks the input directory for a new job (submitted as a json file). 
            Sets the _current_job attribute with the measurement details dictionary.
        """
        while True:
            for file in self._input_dir.glob("*.json"):
                with open(file, "r") as f:

                    job = json.load(f)

                file.unlink()
                self._current_job = job
                return

            time.sleep(5)

    def _prepare_run(self) -> None:
        """
        Simulates the preparation of an HPLC run on the instrument.
        """
        print(f"Sample {self._current_job['name']} loaded.")
        self.silacommunicate(f"Sample {self._current_job['name']} loaded.")

        self._status = "busy"
        print("Preparing the HPLC-MS run")
        self.silacommunicate("Preparing the HPLC-MS run")
        time.sleep(5)
        self._status = "ready"
        print("HPLC-MS is ready for injection.")
        self.silacommunicate("HPLC-MS is ready for injection.")
        

    def _hplc_run(self) -> None:
        """
        Simulates an HPLC-MS run by waiting for injection and waiting for completion of the run.
        """

        time.sleep(5)
        self._status = "busy"
        print("Injection completed. HPLC-MS run started.")

        self.silacommunicate("Injection completed. HPLC-MS run started.")
        time.sleep(10)

        print("HPLC-MS measurement completed.")
        self.silacommunicate("HPLC-MS measurement completed")

        ## 50% chance of crashing
        if random.uniform(0, 1) < 0.4:
            # crash the hplc
            self.client_socket.close()
            sys.exit()


    def _data_analysis(self) -> None:
        print("Data analysis started.")
        self.silacommunicate("Starting data analysis")
        self._current_job["concentration"] = 1.0  # Dummy analysis that always results in a successful identification of the compound.


        if random.uniform(0, 1)< 0.25:
            self.silacommunicate("not detected")
            self._current_job["result"] = "not detected"
        else:
            self.silacommunicate("detected")
            self._current_job["result"] = "detected"

        with open(self._output_dir/f"{self._current_job['name']}.json", "w") as f:
            json.dump(self._current_job, f)

        # TODO: Return an actual data file. 
        shutil.copyfile(os.path.join(STORAGE, "dummy.7z"), os.path.join(OUTPUT, f"{self._current_job['name']}.7z"))


        self.silacommunicate("run completed")

        

        self._current_job = None
        self._status = "idle"
    

    def blank_run(self):
        """
            Simulates a blank 
        """

        self.silacommunicate("run starting")
        time.sleep(5)
        self.silacommunicate("run completed")
        

        self.silacommunicate("Starting data analysis")
        self._current_job["concentration"] = 1.0  # Dummy analysis that always results in a successful identification of the compound.

        with open(self._output_dir / f"{self._current_job['name']}.json", "w") as f:
            json.dump(self._current_job, f)

        # TODO: Return an actual data file. 
        shutil.copyfile(os.path.join(STORAGE, "dummy.7z"), os.path.join(OUTPUT, f"{self._current_job['name']}.7z"))

        self._current_job = None
        self._status = "idle"
    
    
    def silacommunicate(self, message):
        HEADER_LENGTH = 10
        # Prepare username and header and send them
        # We need to encode username to bytes, then count number of bytes and prepare header of fixed size, that we encode to bytes as well
        message = message.encode('utf-8')
        message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
        self.client_socket.send(message_header + message)
        
        


