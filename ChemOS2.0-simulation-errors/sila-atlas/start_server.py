import os

import subprocess
import multiprocessing

import time
from pathlib import Path
import threading
import os
from SilaAtlas import Server


runs = os.path.abspath("atlasruns")
for root, dirs, files in os.walk(runs):
    for f in files:
        os.unlink(os.path.join(root, f))

runs = os.path.abspath("campaigns")
for root, dirs, files in os.walk(runs):
    for f in files:
        os.unlink(os.path.join(root, f))


server_cmd = f"python3 -m SilaAtlas --ip-address 127.0.0.1 -p 65100 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

while True:
    time.sleep(1)

