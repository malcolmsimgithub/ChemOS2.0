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
from utils import run_atlas_calculation, run_chemspeed_synthesis, run_optics_table, run_chemspeed_filter, fetch_hplc_results, fetch_optics_results, hplc_and_chemspeed_inject
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


A001 = 'CC1(C)OB(c2ccc(-n3c4ccccc4c4ccccc43)cc2)OC1(C)C'
B001 = 'CN1CC(=O)OB(OC(=O)C1)\C=C\Br'

engine = create_engine("postgresql://maozer:maozer@localhost:5432/ifogchem")
Session = sessionmaker(bind=engine)
session = Session()



############################
# Create Olympus Campaign  #
############################

# open dict of molecules and ids
with open("job_files/molecules.json", "r") as f:
    molecules = json.load(f)

# create array of fingerprints
fplist = []
for key, value in  molecules.items():
    mol = Chem.MolFromSmiles(value)
    vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048).ToList()
    fplist.append(vector)
# -----------------------------
# Create search space
# -----------------------------
campaign = Campaign()
param_space = ParameterSpace()
param = ParameterCategorical(
    name="building_block", options = list(molecules.keys()), descriptors=fplist
)
param_space.add(param)
campaign.set_param_space(param_space)

campaign.add_observation(["C013"],np.array([0.635E-16, 4.79035141058574e-28]))
campaign.add_observation(["C039"],np.array([0.689E-16, 1.966316933670173e-28]))


with open("laser/0.pkl", "wb") as f:
    pickle.dump(campaign, f)
with open(f"laser/0.pkl", "rb") as f:
        campaign_data = f.read()
    
###########################
# Connect to Sila-atlas   #
###########################

with open("job_files/tanimoto_campaign.json", "r") as f:
    atlas_config = json.load(f)

print("running optimizer")
reccomends = run_atlas_calculation(campaign_data, atlas_config, ATLAS_IP, 65100)

print(reccomends[0])
product = molecules[reccomends[0]]

print(product)

###########################
# Filter crude mixture    #
###########################

instance = run_chemspeed_filter("RACKL:65", CHEMSPEED_IP, 65002, None)
injectposition = instance.get_responses().Termination
print(injectposition)
time.sleep(10)


##########################################################
# run Characterization_1st job on  HPLC                  #
##########################################################

hplc_injection = {
        'name' : f"BSBCz_{reccomends[0]}_jul20", 
        'itype' :'characterization_1st',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : product},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, "SPE_D:3")

while hplc_thread.is_alive():
    time.sleep(1)


#####################################
# Check HPLC results                #
#####################################


job_id = resultsdict["id"]

if fetch_hplc_results(job_id):
    # for file in os.listdir(f"bigjobresults/hplc/{job_id}"):
    #     filename = os.fsdecode(file)
    #     if filename.endswith("results.pkl"):
    #         with open(f"bigjobresults/hplc/{job_id}/{filename}", "rb") as f:
    #             dict = pickle.load(f)

    # if "concentration" in dict.keys():
    #     print("successfully detected compound")
    # else:
    #     print("return compound not detected")
    print("return data found")
else:
    print("return data not found")

###################################################
# Run characterization_2nd Job on HPLC            #
###################################################


hplc_injection = {
        'name' : f"BSBCz_{reccomends[0]}", 
        'itype' :'characterization_2nd',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : product},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, "SPE_D:3")


while hplc_thread.is_alive():
    time.sleep(1)


time.sleep(10)
##########################################################################
# Check HPLC results - get position of compound for optics               #
##########################################################################


job_id = resultsdict["id"]
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
                print(f"Quatum yield: {qy}")
        if f.endswith("AbsPL.csv"):
            with open(f"bigjobresults/optics/{job_id}/{f}", newline='') as file:
                # returning from 2nd row
                data= list(reader(file, delimiter=','))[1]
                tau = float(data[1])
                print(f"tau: {tau}")
        
    
    cross_section_gain = gain*qy/(tau*1e-9*2.25)

print(cross_section_gain)

campaign.add_observation(reccomends[0],np.array([cross_section_gain, 3.144516539068432e-28]))
with open("laser/1.pkl", "wb") as f:
    pickle.dump(campaign, f)