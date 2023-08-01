import subprocess
import time
from pathlib import Path
import threading
import os
from run_socket_server import run_socket_server
from Auto_opt_measurements_ocean import Optical_measurements


#################################
# START SOCKET SERVER           #
#################################
socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()

#################################
# INITIALIZE OPTICS SILA SERVER #
#################################

server_cmd = f"python -m sila2OpticsTable --ip-address 127.0.0.1 -p 65020 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

#################################################
# START OPTICS TABLE MANAGER                    #
#################################################

Opt = Optical_measurements(TE=True, AbsPL= True, Pump= True, Evap=False, logger=True, power_control=True)
Opt.auto_measurement(measurements = ['absorption', 'PL', 'TE'], solvent= 'ACN', redissolution=False, measure_blank=False)