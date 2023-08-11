#!/usr/bin/env python3

from typing import NamedTuple
from pathlib import Path
import json

class Compound(NamedTuple):
    position: int
    volume: float

class PotCfg(NamedTuple):
    v_min: float
    v_max: float
    cycles: int
    steps: int

class ExpCfg(NamedTuple):
    mixing_time: float
    ligand: Compound
    metal: Compound

class RunCfg(NamedTuple):
    potentiostat: PotCfg
    experiment: ExpCfg


def load_cfg_pot(
    dict_cfg: dict
) -> PotCfg:
    return PotCfg(**dict_cfg)


def load_cfg_exp(
    dict_cfg: dict
) -> ExpCfg:
    return ExpCfg (
        mixing_time = dict_cfg["mixing_time"]
        , ligand = Compound(**dict_cfg["ligand"])
        , metal = Compound(**dict_cfg["metal"])
    )

def load_cfg(
    dict_cfg: str
) -> RunCfg:
    return RunCfg (
        potentiostat = load_cfg_pot(dict_cfg["potentiostat"])
        , experiment = load_cfg_exp(dict_cfg["experiment"])
    )

# def run_exp(
#     cfg: ExpCfg
# ) -> None:
#     exp = AutoComplex()
#     exp.run_complexation(
#         num_metal=cfg.metal.position
#         , num_ligand=cfg.ligand.position
#         , quantity_metal=cfg.metal.volume
#         , quantity_ligand=cfg.ligand.volume
#     )

if __name__ == "__main__":
    with open("scratch.json", 'r') as infile:
        json_str = infile.read()
    cfg = load_cfg(json.loads(json_str))
    # run_exp(cfg.experiment)
