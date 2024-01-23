import os

import subprocess
import multiprocessing

import time
from pathlib import Path
import threading
import os
from silasensor import Server


server_cmd = f"python3 -m silasensor --ip-address 127.0.0.1 -p 65100 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

while True:
    time.sleep(1)

