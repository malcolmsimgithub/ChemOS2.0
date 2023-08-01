import subprocess
import threading
from SilaOpticsTable import Server
from run_socket_server import run_socket_server
from optical_characterization.Auto_opt_measurements_ocean import Optical_measurements
import os, shutil


OUTPUT = os.path.abspath("output_folder")

for root, dirs, files in os.walk(OUTPUT):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))


OUTPUT = os.path.abspath("completed_folder")

for root, dirs, files in os.walk(OUTPUT):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))


#################################
# Run Socket server             #
#################################
socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()


#################################
# INITIALIZE Optics Sila Server #
#################################

server_cmd = f"python3 -m SilaOpticsTable --ip-address 127.0.0.1 -p 65070 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)


#################################################
#Boot optics table in simulation mode           #
#################################################

Opt = Optical_measurements(TE=True, AbsPL= True, Pump= True, Evap=False, logger=True, power_control=True)
Opt.auto_measurement(measurements = ['absorption', 'PL', 'TE'], solvent= 'ACN', redissolution=False, measure_blank=False)