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

###########################
# Filter crude mixture    #
###########################


instance = run_chemspeed_filter("RACKL:65", CHEMSPEED_IP, CHEMSPEED_PORT, None)
injectposition = instance.get_responses().Termination
print(injectposition)
time.sleep(10)

###################################################
# Run characterization_2nd Job on HPLC            #
###################################################


hplc_injection = {
        'name' : f"test_injection", 
        'itype' :'characterization_2nd',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : "CCC"},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, injectposition)

while hplc_thread.is_alive():
    time.sleep(1)
print(resultsdict)

if resultsdict["response"].Termination == "HPLCMS_lost":
    while True:
        answer = input("Continue? if so, please restart the hplc and confirm")
        if answer.lower() in ["y","yes"]:
            break
        elif answer.lower() in ["n","no"]:
            continue
        else:
            print("please answer yes or no")
            continue
    hplc_blank()
    print("blank done. re-injecting.")
    hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, injectposition)
    while hplc_thread.is_alive():
        time.sleep(1)
elif resultsdict["response"].Termination ==  "not detected":
    print("compound not detected. rerunning job")
    hplc_thread, chemspeed_thread, resultsdict = hplc_and_chemspeed_inject(hplc_injection, injectposition)
    while hplc_thread.is_alive():
        time.sleep(1)
    print(resultsdict["response"].Termination)
elif resultsdict["response"].Termination ==  "detected":
    print("compound detected!")


