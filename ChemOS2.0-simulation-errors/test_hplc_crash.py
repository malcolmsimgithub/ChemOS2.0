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


hplc_injection = {
        'name' : f"test_injection", 
        'itype' :'characterization_2nd',
        'targets' : [
            {'name' : 'BSBCz_derivative', 'type' : 'product', 'smiles' : "CCC"},
            {'name' : 'BHT', 'type' : 'others', 'formula' : 'C15H24O'}
        ]
}

id_dict = {}

try:
    result = hplc_from_chemspeed(hplc_injection, "SPE_C:9", HPLC_IP, HPLC_PORT, HPLC_CERT, id_dict)
except SilaError:
    while True:
        answer = input("Continue? if so, please restart the hplc and confirm")
        if answer.lower() in ["y","yes"]:
            break
        elif answer.lower() in ["n","no"]:
            sys.exit()
        else:
            print("please answer yes or no")
            continue
    try:
        print("performing blank run on HPLCMS...")
        hplc_blank()
        result = hplc_from_chemspeed(hplc_injection, "SPE_C:9", HPLC_IP, HPLC_PORT, HPLC_CERT, id_dict)
    except SilaError:
        print("HPLC has crashed twice. please inspect for detailed measurements")
        sys.exit()

print(result)

if result ==  "not detected":
    print("compound not detected. rerunning job")
    try:
        result = hplc_from_chemspeed(hplc_injection, "SPE_C:9", HPLC_IP, HPLC_PORT, HPLC_CERT, id_dict)
    except SilaError:
        print("HPLC has crashed. please inspect for detailed measurements or run again")
        sys.exit()
elif result ==  "detected":
    print("compound detected!")

sys.exit()