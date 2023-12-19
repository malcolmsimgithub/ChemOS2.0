
import subprocess
from chemspeed_operator_process.Operator import ChemSpeedOperator
from pathlib import Path
import threading
from run_socket_server import run_socket_server
import os, shutil

#######################################
# SIMULATION MODE - SET THIS VARIABLE #
#######################################

SIMULATION = True

############################################################
# DEFINITION OF PATHS TO INSTRUMENTS AND IMPORTANT FOLDERS #
############################################################

if SIMULATION:
    CONTROLLER: Path = Path(__file__).parent
else:
    CONTROLLER: str = r'\\IPC5\Users\Operator\Desktop\Commands'

SETTINGS: Path = Path("chemspeed_operator_process/Default_Settings")
OUTPUT: Path = Path("chmspd_output")

CHAR = OUTPUT/"Characterization/Completed"
SYNTH = OUTPUT/"Synthesis/Completed_Batches"
TOMAKE = OUTPUT/"Synthesis/Batches_to_Make"

for root, dirs, files in os.walk(TOMAKE):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))



for root, dirs, files in os.walk(CHAR):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))

for root, dirs, files in os.walk(SYNTH):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))

SYNTHESIS_PROCEDURE: Path = OUTPUT/"Synthesis/Two_Step_Suzuki_TIDA.json"


#################################
# INITIALIZE CHEMSPEED OPERATOR #
#################################

server_cmd = f"python -m chmspd_sila2_pkg --ip-address 127.0.0.1 -p 65001 --insecure"

proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)


socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()

chemspeed = ChemSpeedOperator(command_path=CONTROLLER,
                              defaults_folder=SETTINGS,
                              output_path=OUTPUT,
                              clear_folders=True,
                              threads=("Synthesis", "Characterization"),
                              simulation=SIMULATION,
                              auto_start=True
                              )
