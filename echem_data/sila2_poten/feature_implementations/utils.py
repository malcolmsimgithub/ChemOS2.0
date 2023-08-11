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
    num_mixings: int
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
        num_mixings = dict_cfg["num_mixings"]
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



from pathlib import Path
import numpy as np
from numpy.typing import NDArray
import pandas as pd
import sys
from scipy.signal import find_peaks

def process_cycles(
    csv_path: Path
):
    global reduction, oxidation
    HEADER = ("I", "V", "R", "T", "h", "t", "c")
    df = pd.read_csv(
        csv_path
        , names=HEADER
        , index_col=False
    )
    df = df[-0.75 < df["V"]]
    df = df[df["V"] < 0.75]

    df["c"] = df["c"].astype(int)
    # Get Red Ox
    (_, reduction), (_, oxidation) = df.groupby(np.diff(df["V"].array, append=0) > 0)
    oxi_peaks = map(
        lambda xs: xs[1][["V", "I"]].iloc[find_peaks(xs[1]["I"], width=10, distance=10, rel_height=0.9)[0]]
        , oxidation.groupby("c")
    )
    red_peaks = map(
        lambda xs: xs[1][["V", "I"]].iloc[find_peaks(-xs[1]["I"], width=10, distance=10, rel_height=0.2)[0]]
        , reduction.groupby("c")
    )
    return tuple(oxi_peaks), tuple(red_peaks)


if __name__ == "__main__":
    oxi_peaks, _ = process_cycles(sys.argv[1])
    def checknum(x):
        y = x["I"].argmax()
        if (pd.isna(y)): return -1
        return x["V"].iloc[y]

    try:
        oxivalues = list(map(checknum,oxi_peaks))[1:]
        max_val= np.max(oxivalues)

    except:
        max_val = -5

    print(max_val)


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
