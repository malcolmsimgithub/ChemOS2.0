import os
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from olympus.objects import ParameterCategorical
import pathlib
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

BUDGET=3
PATH: pathlib.Path = pathlib.Path().resolve()
RUNSDIR = os.path.join(PATH, "examples/dynamic_search_space/runs")
MODELS = [
    "RandomSearch",
    "Botorch",
    "RGPE",
    "DKT",
    "Dynamic"
]

# open dict of molecules and ids
with open("molecules.json", "r") as f:
    molecules = json.load(f)


# create array of fingerprints
fplist = []
for key, value in  molecules.items():
    mol = Chem.MolFromSmiles(value)
    vector = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048).ToList()
    fplist.append(vector)


# -----------------------------
# Instantiate surface
# -----------------------------


campaign = Campaign()
param_space = ParameterSpace()
# add 3 continuous Parameters
param = ParameterCategorical(
    name="building_block", options = list(molecules.keys()), descriptors=fplist
)

param_space.add(param)
campaign.set_param_space(param_space)

with open("sila.pkl", "wb") as f:
    pickle.dump(campaign, f)