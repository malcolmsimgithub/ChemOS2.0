from sila2.client import SilaClient
import time
from pathlib import Path
from sila2.client import SilaClient
import time
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from olympus.objects import ParameterCategorical
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import json
from database import *
from sqlalchemy import (
    create_engine,
)
from sqlalchemy.orm import sessionmaker
import json
from utils import *
import pickle
from csv import reader
import sys

###########################
# Begin Characterization  #
###########################
def timestamp():
    return time.strftime("%y-%m-%d-%H-%M-%SS", time.localtime())

with open("job_files/optical_job.json", "r") as f:
    opticsjob = json.load(f)

opticsjob["name"] = f"BSBCz_derivative_{timestamp()}"

try:
    message, job_id =  run_optics_table(opticsjob, OPTICS_IP, OPTICS_PORT, OPTICS_CERT)
except SilaError:
    while True:
        answer = input("Optics table has crashed! if so, please restart the optics table and confirm")
        if answer.lower() in ["y","yes"]:
            break
        elif answer.lower() in ["n","no"]:
            continue
        else:
            print("please answer yes or no")
            continue
    print("re-doing optics table job")
    opticsjob["name"] = f"BSBCz_derivative_{timestamp()}"
    try:
        message, job_id =  run_optics_table(opticsjob, OPTICS_IP, OPTICS_PORT, OPTICS_CERT)
    except SilaError:
        print("optics table has failed twice. aborting")
        sys.exit()

print("waitng 10 seconds before analysis of results")
time.sleep(10)

resultsfound = fetch_optics_results(job_id)
if resultsfound:
    for f in os.listdir(f"jobresults/optics/{job_id}"):
        if f.endswith("AbsPL.csv"):
            with open(f"jobresults/optics/{job_id}/{f}", newline='') as file:
                # returning from 2nd row
                data= list(reader(file, delimiter=','))[1]
                gain = float(data[11])
                print(f"optical gain factor: {gain}")

                qy = float(data[10])
                print(f"Quantum yield: {qy}")
        if f.endswith("TE.csv"):
            with open(f"jobresults/optics/{job_id}/{f}", newline='') as file:
                # returning from 2nd row
                data= list(reader(file, delimiter=','))[1]
                tau = float(data[1])
                print(f"tau: {tau}")
        
    
    cross_section_gain = gain*qy/(tau*1e-9*2.25)
print(cross_section_gain)
