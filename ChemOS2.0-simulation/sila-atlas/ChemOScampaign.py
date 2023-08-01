import os
import numpy as np
from olympus import Campaign
from olympus.campaigns import ParameterSpace
from atlas.optimizers.base import BasePlanner
from atlas.optimizers.gp import BoTorchPlanner
from atlas.optimizers.tanimoto import TanimotoPlanner
# from atlas.optimizers.dynamic_space.planner import DynamicSSPlanner
from atlas.optimizers.params import Parameters, ParameterContinuous, ParameterCategorical
import pathlib
import json
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

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


# for x in X_init:
#     measurement = surface(x.reshape((1, x.shape[0])))
#     campaign.add_observation(x, measurement)

value_space = ParameterSpace()
value_space.add(ParameterContinuous(name='Simulation gain'))
value_space.add(ParameterContinuous(name='Experimental gain'))



planner = TanimotoPlanner(
    ############ TODO DESIGN STRATEGY
    goal='minimize',
    init_design_strategy='random',
    num_init_design=1,
    batch_size=1,
    use_descriptors=True,
    acquisition_optimizer_kind='genetic',
    is_moo=True,
    scalarizer_kind='Hypervolume', 
    ###### TODO VALUE SPACE
    value_space=value_space,
    goals=['max', 'max'],   
)

planner.set_param_space(param_space)



# start the optimization experiment
iteration = 0
# optimization loop
while len(campaign.values) < BUDGET:
    print("############################")

    print(f"\nITERATION : {iteration+1}\n")

    samples = planner.recommend(campaign.observations)
    print(f"SAMPLES : {samples}")

    for sample in samples:
        sample_arr = sample.to_array()
        # measurement = surface.run(
        #     sample_arr.reshape((1, sample_arr.shape[0]))
        # )
        
        campaign.add_observation(sample_arr, np.array([15,63]))

        print(campaign.params)
        print(campaign.values)

    
    iteration += 1
    print("############################")


print(f"run completed")
