from aiida.engine import run
from aiida.orm import Str
from laser_workchain import SilaLaserWorkChain

with open("input.txt", "r") as f:
    smiles = f.read()
run(SilaLaserWorkChain, smiles=Str(smiles))