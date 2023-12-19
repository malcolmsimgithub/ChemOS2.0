import pickle

import subprocess

with open("simulator_storage/20201204/1_c_1_vial5.pkl", "rb") as f:
    data = pickle.load(f)


cmd = f"7z a test/20201204.7z simulator_storage/20201204"
proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)