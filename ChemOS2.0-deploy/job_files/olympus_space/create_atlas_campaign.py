import os
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from olympus.objects import ParameterCategorical, ParameterContinuous
import pathlib
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

BUDGET=3
PATH: pathlib.Path = pathlib.Path().resolve()

# -----------------------------
# Instantiate surface
# -----------------------------


campaign = Campaign()
param_space = ParameterSpace()
# add 3 continuous Parameters
param = ParameterCategorical(
    name="ligand_identity", options = ["","",""], descriptors=None
)
param_space.add(param)

param = ParameterContinuous(name="mixing_time", high=60, low=0
)
param_space.add(param)

param = ParameterContinuous(name="ligand/metal ratio", high=0.99, low=0.01
)
param_space.add(param)
campaign.set_param_space(param_space)

with open("sila.pkl", "wb") as f:
    pickle.dump(campaign, f)