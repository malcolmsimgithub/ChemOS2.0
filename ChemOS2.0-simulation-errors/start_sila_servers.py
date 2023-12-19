import subprocess
import time


server_cmd = f"cd sila-atlas; python3 start_server.py"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)


server_cmd = f"cd sila-optics; python3 start_server.py"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

server_cmd = f"cd sila-hplc; python3 start_server.py"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

server_cmd = f"cd sila-chemspeed; python3 start_server.py"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

while True:
    time.sleep(1)