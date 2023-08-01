import time
import psycopg2 as psycopg
from pathlib import Path
import os
import time
import json
from connection import *
from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    ForeignKey,
    DateTime,
    LargeBinary
)
from sqlalchemy.orm import sessionmaker
import json
import pickle
from csv import reader
from utils import run_atlas_calculation, run_chemspeed_synthesis, run_optics_table, run_chemspeed_filter, fetch_hplc_results, fetch_optics_results, hplc_from_chemspeed, hplc_from_autosampler
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from olympus.objects import ParameterCategorical
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
import numpy as np

CHEMSPEED_IP = "127.0.0.1"
HPLC_IP = "127.0.0.1"
OPTICS_IP = "127.0.0.1"
ATLAS_IP = "127.0.0.1"


smiles = "CC1(C)c2cc(/C=C/c3ccc(-n4c5ccccc5c5ccccc54)cc3)ccc2-c2ccc(/C=C/c3ccc(-n4c5ccccc5c5ccccc54)cc3)cc21"

engine = create_engine("postgresql://malcomr:malcolm@localhost:5432/ifogchem")
Session = sessionmaker(bind=engine)
session = Session()


##########################################################
# run Characterization_1st job on  HPLC                  #
##########################################################

hplc_injection = {
        'name' : f"BSBCz_C005", 
        'itype' :'characterization_1st',
        'position': "Y:A2",
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : smiles},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

job_id = hplc_from_autosampler(hplc_injection, HPLC_IP, 65004, None)


#####################################
# Check HPLC results                #
#####################################

if fetch_hplc_results(job_id):
    for file in os.listdir(f"bigjobresults/hplc/{job_id}"):
        filename = os.fsdecode(file)
        if filename.endswith("results.pkl"):
            with open(f"bigjobresults/hplc/{job_id}/{filename}", "rb") as f:
                dict = pickle.load(f)
else:
    print("return data not found")

###################################################
# Run characterization_2nd Job on HPLC            #
###################################################


hplc_injection = {
        'name' : f"BSBCz_C005", 
        'position': "Y:A2",
        'itype' :'characterization_2nd',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : smiles},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}


job_id = hplc_from_autosampler(hplc_injection, HPLC_IP, 65004, None)



##########################################################################
# Check HPLC results - get position of compound for optics               #
##########################################################################

if fetch_hplc_results(job_id):
    with open(f"bigjobresults/hplc/{job_id}/characterization_params.json", "r") as f:
        opticsjob = json.load(f)



###########################
# Begin Characterization  #
###########################
def timestamp():
    return time.strftime("%y-%m-%d-%H-%M", time.localtime())


opticsjob["name"] = f"BSBCz_derivative_{timestamp()}"

instance, job_id =  run_optics_table(opticsjob, OPTICS_IP, 65017, None)

while instance.done != True:
    time.sleep(1)

print("waitng 10 seconds before analysis of results")
time.sleep(10)

resultsfound = fetch_optics_results(job_id)

if resultsfound:
    for f in os.listdir(f"bigjobresults/optics/{job_id}"):
        if f.endswith("AbsPL.csv"):
            with open(f"bigjobresults/optics/{job_id}/{f}", newline='') as file:
                # returning from 2nd row
                data= list(reader(file, delimiter=','))[1]
                gain = float(data[11])
                print(f"optical gain factor: {gain}")

                qy = float(data[10])
                print(f"Quantum yield: {qy}")
        if f.endswith("TE.csv"):
            with open(f"bigjobresults/optics/{job_id}/{f}", newline='') as file:
                # returning from 2nd row
                data= list(reader(file, delimiter=','))[1]
                tau = float(data[1])
                print(f"tau: {tau}")
        
   
