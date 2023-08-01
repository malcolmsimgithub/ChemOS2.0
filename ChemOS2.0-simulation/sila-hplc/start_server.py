from hplcsimulator import HPLCSimulator
from pathlib import Path
import subprocess
import threading
from socket_server import run_socket_server
import os
import shutil
COMMDIR :Path = Path("communication/")
OUTPUT = os.path.join(COMMDIR, "jobs_returned")


for root, dirs, files in os.walk(OUTPUT):
    for f in files:
        os.unlink(os.path.join(root, f))
    for d in dirs:
        shutil.rmtree(os.path.join(root, d))


server_cmd = f"python3 -m silahplc --ip-address 127.0.0.1 -p 65010 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()

simulator = HPLCSimulator(COMMDIR)
simulator.operate()







