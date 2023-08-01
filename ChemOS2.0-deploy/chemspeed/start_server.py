import os
import logging
import subprocess
import multiprocessing
from chemspeed_operator_process.Operator import ChemSpeedOperator
import time
from pathlib import Path
import threading
import os
from chmspd_sila2_pkg import Server
from run_socket_server import run_socket_server


# SIMULATION MODE - SET THIS VARIABLE #
#######################################

SIMULATION =True

############################################################
# DEFINITION OF PATHS TO INSTRUMENTS AND IMPORTANT FOLDERS #
############################################################

if SIMULATION:
    CONTROLLER: Path = Path(__file__).parent
else:
    CONTROLLER: str = r'\\IPC5\Users\Operator\Desktop\Commands'

#DROPBOX_PATH = get_dropbox_path() / "PythonScript"
SETTINGS: Path = "chemspeed_operator_process/Default_Settings"
OUTPUT: Path = "chmspd_output"
SYNTHESIS_PROCEDURE: Path = os.path.join(OUTPUT,"Synthesis/Two_Step_Suzuki_TIDA.json")



#######################
# START SILA2 SERVER  #
#######################
server_cmd = f"python -m chmspd_sila2_pkg --ip-address 127.0.0.1 -p 65040 --cert-file chemspeed.crt -k chemspeed.key"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)


#################################
# INITIALIZE SOCKET SERVER      #
#################################
socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()


#################################
# INITIALIZE CHEMSPEED OPERATOR #
#################################

chemspeed = ChemSpeedOperator(command_path=CONTROLLER,
                              defaults_folder=SETTINGS,
                              output_path=OUTPUT,
                              clear_folders=False,
                              threads=("Synthesis", "Characterization"),
                              #synthesis_procedure=SYNTHESIS_PROCEDURE,
                              simulation=SIMULATION,
                              auto_start=True
                              )
