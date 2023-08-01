import threading
import time
from pathlib import Path
from typing import Union
from singleton_decorator import singleton
import socket

from .Executor import ChemSpeedExecutor
from .ChemSpeedModules import ChemSpeedFactory
from .Workflows import Batch, CompoundStorage, HPLCCharacterizationScheduler
from .Utils import UserInterface, load_json
import os


IP = "127.0.0.1"
PORT = 65002
HEADER_LENGTH = 10
OUTPUT: Path = "chmspd_output"
PROCEDUREDIR = os.path.join(os.path.join(OUTPUT, "Synthesis"), "Procedure")

def status_active():
        statusfile = os.path.join(OUTPUT, "status.txt")
        with open(statusfile, "w") as f:
            f.write("active")
def status_inactive():
    statusfile = os.path.join(OUTPUT, "status.txt")
    with open(statusfile, "w") as f:
        f.write("inactive")


@singleton
class ChemSpeedOperator(object):
    """
    Global Operator Object for Complex Synthesis and Characterization Workflows using the ChemSpeed Platform.
        -   Instrument configuration is specified via the HardwareConfiguration.json file.
        -   Manages batch-wise synthesis (via a general json synthesis procedure and batch specifications) and
            characterization via the HPLC injection port.
        -   Opens a TK user interface (mostly for visualization / logging) via multithreading.

    Public methods:
        start_operation(operations: tuple, synthesis_procedure: Path):  Starts the operation and opens the TK
                                                                        user interface
                                                                        Is automatically called when auto_start=True
    """

    def __init__(
            self,
            command_path: Path,
            defaults_folder: Path,
            output_path: Path,
            clear_folders: bool,
            threads: Union[tuple, str],
            synthesis_procedure: Path = None,
            simulation: bool = False,
            auto_start: bool = True,
            **kwargs):
        """
        Instantiates the ChemSpeed Operator object.

        Args:
            command_path: Path to the ChemSpeed Folder
            defaults_folder: Path to the folder where the default settings are located.
            output_path: Path to the folder where outgoing settings are located.
            clear_folders: Whether to clear existing settings (compound storage, filtration positions, ...) from folders
            threads: Operations to be performed: ("Synthesis", "Characterization") or only one of these options.
            synthesis_procedure: Path to the file where the synthesis procedure .json is located
            simulation: Simulation mode (True or False)
            auto_start: Start all requested threads upon instantiation (True or False)
            kwargs: Keyword arguments for the user interface to use.
        """
        status_inactive()
        self.paths = self._setup_folders(output_path, defaults_folder, clear_folders)

        self.hardware = dict()
        self._initialize_hardware()

        self.silasocket = client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Connect to a given ip and port
        self.silasocket.connect((IP, PORT))

        # Set connection to non-blocking state, so .recv() call won;t block, just return some exception we'll handle
        self.silasocket.setblocking(False)

        id  = 'chemspeed'.encode('utf-8')
        username_header = f"{len(id):<{HEADER_LENGTH}}".encode('utf-8')
        client_socket.send(username_header + id)

        self.user_interface = UserInterface(self.silasocket,threads=threads,**kwargs)

        self.storage = CompoundStorage(
            components=self.hardware,
            output_path=self.paths["storage"],
            defaults_path=self.paths["defaults"],
            ui_confirmation=self.user_interface.request_confirmation
        )

        self.executor = ChemSpeedExecutor(
            command_path=command_path,
            output_path=self.paths["output"],
            compound_location_mapping=self.storage.get_compound_location,
            reactor=self.hardware["reactor"],
            simulation=simulation
        )


        self.simulation = simulation

        if auto_start:
            self.start_operation(threads, synthesis_procedure)

    @staticmethod
    def _setup_folders(
            parent_dir: Path,
            defaults_dir: Path,
            clear_folders: bool
    ) -> dict:
        """
        Sets up the (output) folder architecture required for the ChemSpeedOperator, if non-existent.
        Depends on the parent directory:

            /Synthesis/
                /Batches_to_Make/
                /Completed_Batches/
            /Characterization/
            /ChemSpeed_Modules/
            /Storage_Inventory/

        Args:
            parent_dir: Parent Directory
            defaults_dir: Directory of default files
            clear_folders: Whether to clean up the folders "modules", "storage" and "characterization"

        Returns:
            paths: Dictionary of path names and corresponding Path objects.
        """
        paths = {
            "defaults": defaults_dir,
            "output": parent_dir,
            "synthesis": parent_dir/"Synthesis",
            "batches_to_make": parent_dir/"Synthesis/Batches_to_Make",
            "completed_batches": parent_dir/"Synthesis/Completed_Batches",
            "characterization":parent_dir/"Characterization/Characterizations_to_do",
            "characterization_output":parent_dir/"Characterization/Completed",
            "modules": parent_dir/"ChemSpeed_Modules",
            "storage": parent_dir/"Storage_Inventory",
        }

        for path in paths.values():
            if os.path.isdir(path)==False:
                os.mkdir(path)

        if clear_folders:
            for folder in ("characterization", "modules", "storage"):
                for file in paths[folder].iterdir():
                    if file.is_file():
                        file.unlink()

        return paths

    def _initialize_hardware(
            self
    ) -> None:
        """
        Initialize the self.hardware attribute from the corresponding settings file.
        Generates the dictionary with all instantiated hardware objects.
        """
        hardware_settings = load_json(os.path.join(self.paths["defaults"], "HardwareConfig_old.json"))
        # hardware_settings = load_json(os.path.join(self.paths["defaults"], "HardwareConfig.json"))
        hardware_factory = ChemSpeedFactory()

        for key, module_settings in zip(hardware_settings, hardware_settings.values()):
            self.hardware[key] = hardware_factory(
                object_type=module_settings["type"],
                output_path=self.paths["output"],
                defaults_path=self.paths["defaults"],
                **module_settings
            )

    def start_operation(
            self,
            operations: Union[str, tuple],
            synthesis_procedure: Path
    ) -> None:
        """Starts the operation, as specified in the threads list ("Synthesis", "Characterization").
        Creates a thread for each operation.
        Starts the graphical user interface.
        ATTENTION: Blocks the main thread

        Args:
            operations: Workflows to be run: ("synthesis", "characterization") or one of these options as tuple.
            synthesis_procedure: Path to the synthesis procedure to be executed in the self._run_synthesis function

        Returns:
            None
        """
        status_inactive()

        if "Synthesis" in operations:
            threading.Thread(target=self._run_synthesis, kwargs={"synthesis_procedure": synthesis_procedure}).start()
        if "Characterization" in operations:
            threading.Thread(target=self._run_characterization).start()

        threading.Thread(target=self._human_intervention).start()
        self.user_interface.activate_gui()

    #########################################################################
    # HIGH-LEVEL METHODS FOR ENTIRE WORKFLOWS (SYNTHESIS, CHARACTERIZATION) #
    #########################################################################

    def _run_synthesis(
            self,
            synthesis_procedure: Path
    ) -> None:
        """General function that runs a specific synthesis whenever a batch is available in the respective folder.

        Args:
            synthesis_procedure (pathlib.Path): Path to the json file containing the synthesis procedure
        """
        self.user_interface.communicate("Synthesis Channel Started.", thread="Synthesis")

        while True:
            while synthesis_procedure == None:
                if len(os.listdir(PROCEDUREDIR)) == 0:
                    self.executor.status = "interrupt"
                    self.user_interface.communicate(f"Procedure file required before starting synthesis channel", thread="Synthesis")
                    time.sleep(1)
                else: 
                    self.synthesis_procedure = os.path.join(PROCEDUREDIR, os.listdir(PROCEDUREDIR)[0])
                    synthesis_procedure = os.path.join(PROCEDUREDIR, os.listdir(PROCEDUREDIR)[0])
                    old_procedure = synthesis_procedure
                    self.executor.status = "idle"
                    self.user_interface.communicate(f"Procedure file being used is {os.listdir(PROCEDUREDIR)[0]}", thread="Synthesis")
                    break
            self.synthesis_procedure = os.path.join(PROCEDUREDIR, os.listdir(PROCEDUREDIR)[0])
            synthesis_procedure = os.path.join(PROCEDUREDIR, os.listdir(PROCEDUREDIR)[0])
            if synthesis_procedure != old_procedure:
                self.user_interface.communicate(f"Procedure file being used is {os.listdir(PROCEDUREDIR)[0]}", thread="Synthesis")
                old_procedure = synthesis_procedure
            self._synthesize_batch(synthesis_procedure)
            time.sleep(5)

    def _synthesize_batch(
            self,
            procedure: Path
    ) -> None:
        """
        Performs the synthesis of a batch by:
            - Loading the batch from the specified folder (if a batch file is available)
            - Instantiating a Batch object, which modifies the synthesis procedure
            - Parsing the batch by assigning it to specific reaction vials and generating a list of Operations
            - Updating the compound storage object (and the storage racks) accordingly.
            - Executing all required Operations sequentially.
            - Clearing the batch and the storage eventually.

        Args:
            procedure (Path): Path to the synthesis procedure .json file.
        """
        status_active()
        batch = Batch(synthesis_procedure=procedure, synthesis_folder=self.paths["synthesis"])
        if not batch.batch:  # Returns if no batch was loaded
            return

        self.user_interface.communicate(f"Synthesis of batch {batch.batch_name} started.", thread="Synthesis")

        batch.process_batch(reactor=self.hardware["reactor"], step_settings=self.paths["defaults"])
        self.storage.generate_inventories(batch.get_operations_list())

        try:
            while True:
                operation_str, target, source, kwargs = batch.run_next_operation()
                self._run_operation(operation_str, target, source, **kwargs)
        except StopIteration:
            batch.clear_batch()
            self.hardware["reactor"].clear_reactor()
            self.storage.clear_storage()
        status_inactive

    def _run_characterization(
            self
    ) -> None:
        """
        Starts continuous characterization runs (should be run on a single thread).
        Initializes a HPLCCharacterizationScheduler object which generates the operations to be executed.
        """
        characterization_scheduler = HPLCCharacterizationScheduler(
            default_settings=self.paths["defaults"],
            characterization_path=self.paths["characterization"],
            output_path=self.paths["characterization_output"],
            hardware=self.hardware
        )

        self.user_interface.communicate("Characterization Channel Started.", thread="Characterization")

        while True:
            time.sleep(5)
            try:
                operation_str, target_zone, source_zone, kwargs = characterization_scheduler.execute_next_operation()
                self._run_operation(operation_str, target_zone, source_zone, **kwargs)
            except StopIteration:
                continue

    def _human_intervention(
            self
    ) -> None:
        """
        Checks for human intervention via the buttons in the GUI.
        Initiates the corresponding actions:
            - "interrupt" -> instrument set to busy
            - "emergency" -> instrument set to busy, all tools unmounted
            - "shutdown"  -> instrument set to busy, all tools unmounted, manger shut down
            - "idle"      -> instrument status set back to idle, threads should resume
        """
        while True:
            if self.user_interface.gui.status:
                if self.user_interface.gui.status == "interrupt":
                    self.user_interface.communicate("Operation Interrupted!", thread="all")
                    self.executor.status = "interrupt"
                    while self.user_interface.gui.status == "interrupt":
                        time.sleep(1)
                elif self.user_interface.gui.status == "emergency":
                    self.user_interface.communicate("Operation Interrupted! Unmount of all Tools Initiated!", thread="all")
                    self.executor.unmount_all()
                    self.executor.status = "interrupt"
                    while self.user_interface.gui.status == "emergency":
                        time.sleep(1)
                elif self.user_interface.gui.status == "restart":
                    self.user_interface.communicate("Operation Resumed!", thread="all")
                    self.executor.status = "idle"
                    self.user_interface.gui.set_status(None)
                elif self.user_interface.gui.status == "shutdown":
                    self.user_interface.communicate("Operation Interrupted! Unmount of all Tools Initiated! Manager will be shut down!", thread="all")
                    self.executor.unmount_all()
                    self.executor.shutdown_manager()
                    self.executor.status = "interrupt"
                    while self.user_interface.gui.status == "shutdown":
                        time.sleep(1)

    ############################################
    # LOW-LEVEL EXECUTION OF SINGLE OPERATIONS #
    ############################################

    def _run_operation(
            self,
            operation: str,
            target_zone: Union[str, list],
            source_zone: Union[str, None],
            **kwargs
    ) -> None:
        """Parse operation strings from synthesis procedure and call respective methods.

        Args:
            operation: Operation to be performed (needs to be in available_operations)
            target_zone: name of the target_zone zone(s), can be empty
            source_zone: name of the source_zone zone, can be None
            **kwargs:  Keyword Arguments for the Operation to be performed
        """
        available_operations = {
            "transfer_compound": self.executor.transfer_compound,
            "schlenk_cycle": self.executor.schlenk_cycling,
            "reflux": self.executor.reflux,
            "prime_pumps": self.executor.prime_pumps,
            "set_drawer": self.executor.set_drawer,
            "communicate": self.user_interface.communicate,
            "request_confirmation": self.user_interface.request_confirmation,
            "filter_collect": self.executor.filter_collect,
            "inject_to_hplc": self.executor.inject_to_hplc,
        }

        return available_operations[operation](target_zone=target_zone, source_zone=source_zone, **kwargs)

