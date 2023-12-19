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


CHEMSPEED_IP = "127.0.0.1"
CHEMSPEED_CERT = None
CHEMSPEED_PORT = 65001
HPLC_IP = "127.0.0.1"
HPLC_CERT = None
HPLC_PORT = 65010
OPTICS_IP = "127.0.0.1"
OPTICS_CERT = None
OPTICS_PORT = 65070
ATLAS_IP = "127.0.0.1"
ATLAS_PORT = 65100


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


campaign.add_observation(["C039"],np.array([0.689E-16]))
campaign.add_observation(["C013"],np.array([0.635E-16]))


with open("runs/0.pkl", "wb") as f:
    pickle.dump(campaign, f)


for i in range(10):
    with open(f"runs/{i}.pkl", "rb") as f:
        campaign = pickle.load(f)
    with open(f"runs/{i}.pkl", "rb") as f:
        campaign_data = f.read()

    ###########################
    # Connect to Sila-atlas   #
    ###########################

    print("running optimizer")
    reccomends = run_atlas_calculation(campaign_data, atlas_config, ATLAS_IP, 65100)
    product = molecules[reccomends[0]]


    ###########################
    # Create batch file       #
    ###########################

    compounddict = {"Atlassuggestion": {"$A$": "A001", "$B$": "A001", "$C$": reccomends[0], "vial": "ISYNTH:30"}}

    with open('job_files/chemspeed_synthesis.json', 'w') as fp:
        json.dump(compounddict, fp)
        print("successfully optimized and created chemspeed synthesis file")

    ###########################
    # Connect to Chemspeed    #
    ###########################

    file = Path("job_files/chemspeed_synthesis.json")
    f = open(file)
    jsondata = dict(json.load(f))
    instance = run_chemspeed_synthesis(jsondata, "test synthesis", CHEMSPEED_IP, CHEMSPEED_PORT, CHEMSPEED_CERT)

    while instance.done != True:
        time.sleep(1)

    ###########################
    # Filter crude mixture    #
    ###########################


    instance = run_chemspeed_filter("RACKL:65", CHEMSPEED_IP, CHEMSPEED_PORT, None)
    injectposition = instance.get_responses().Termination
    time.sleep(10)

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

    hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, injectposition)


    while hplc_thread.is_alive():
        time.sleep(1)


    time.sleep(10)

    ##########################################################################
    # Check HPLC results - get position of compound for optics               #
    ##########################################################################


    job_id = resultsdict["id"]

    print(f"job id: {job_id}")
    if fetch_hplc_results(job_id):
        with open(f"jobresults/hplc/{job_id}/characterization_params.json", "r") as f:
            opticsjob = json.load(f)

    ###########################
    # Begin Characterization  #
    ###########################
    def timestamp():
        return time.strftime("%y-%m-%d-%H-%M", time.localtime())

    opticsjob["name"] = f"BSBCz_derivative_{timestamp()}"


    instance, job_id =  run_optics_table(opticsjob, OPTICS_IP, OPTICS_PORT, OPTICS_CERT)

    while instance.done != True:
        time.sleep(1)

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

    campaign.add_observation([reccomends[0]],np.array([cross_section_gain]))
    with open(f"runs/{i+1}.pkl", "wb") as f:
        pickle.dump(campaign, f)