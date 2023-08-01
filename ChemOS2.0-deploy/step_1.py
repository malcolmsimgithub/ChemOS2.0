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
from utils import run_atlas_calculation, run_chemspeed_synthesis
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from olympus.objects import ParameterCategorical
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
import numpy as np
import subprocess


CHEMSPEED_IP = "127.0.0.1"
HPLC_IP = "127.0.0.1"
OPTICS_IP = "127.0.0.1"
ATLAS_IP = "127.0.0.1"

A001 = 'CC1(C)OB(c2ccc(-n3c4ccccc4c4ccccc43)cc2)OC1(C)C'
B001 = 'CN1CC(=O)OB(OC(=O)C1)\C=C\Br'

engine = create_engine("postgresql://maozer:maozer@localhost:5432/ifogchem")
Session = sessionmaker(bind=engine)
session = Session()

Aiida_input = '/home/malcolm/sila-aiida/input.txt'
virtual_gain_factor = '/home/malcolm/sila-aiida/final_step/gain_factor.json'


############################
# Create Olympus Campaign  #
############################

# open dict of molecules and ids
with open("job_files/molecules.json", "r") as f:
    molecules = json.load(f)

with open(f"job_files/laser_campaign.json", "r") as f:
        atlas_config = json.load(f)

# create array of fingerprints
fplist = []
for key, value in  molecules.items():
    mol = Chem.MolFromSmiles(value)
    vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048).ToList()
    fplist.append(vector)


# Create search space
campaign = Campaign()
param_space = ParameterSpace()
param = ParameterCategorical(
    name="building_block", options = list(molecules.keys()), descriptors=fplist
)
param_space.add(param)
campaign.set_param_space(param_space)

campaign.add_observation(["C039"],np.array([0.689E-16, 1.966316933670173e-28]))
campaign.add_observation(["C013"],np.array([0.635E-16, 4.79035141058574e-28]))


with open("runs/laser_0.pkl", "wb") as f:
    pickle.dump(campaign, f)


###########################
# Connect to Sila-atlas   #
###########################

with open("runs/laser_0.pkl", "rb") as f:
    campaign = pickle.load(f)
with open("runs/laser_0.pkl", "rb") as f:
    campaign_data = f.read()

print("running optimizer")
reccomends = run_atlas_calculation(campaign_data, atlas_config, ATLAS_IP, 65100)
product = molecules[reccomends[0]]

###########################
# Create batch file       #
###########################

compounddict = {"BSBCz_derivative": {"$A$": "A001", "$B$": "B001", "$C$": reccomends[0], "MW": 598.26}}

with open('job_files/BSBCz_derivative.json', 'w') as fp:
    json.dump(compounddict, fp)
    print("successfully optimized and created chemspeed synthesis file")


with open(Aiida_input, "w") as f:
    f.write(product)

###########################
# Begin DFT workchain     #
###########################

cmd = f"""
            cd /home/malcolm/sila-aiida && verdi run run_workchain.py
            """

proc = subprocess.Popen(cmd, shell=True)


###########################
# Connect to Chemspeed    #
###########################

file = Path("job_files/BSBCz_derivative.json")
f = open(file)
jsondata = dict(json.load(f))
instance = run_chemspeed_synthesis(jsondata, "BSBCz_derivative_synthesis", CHEMSPEED_IP, 65002, None)

while instance.done != True:
    time.sleep(1)


proc.communicate()

with open(virtual_gain_factor, "r") as f:
    file = json.load(f)
    print("virtual gain factor:")
    print(file["max_gain_factor"])