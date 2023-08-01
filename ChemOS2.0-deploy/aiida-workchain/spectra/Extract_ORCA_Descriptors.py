
#!/usr/bin/env python
from typing import Dict, List, Tuple, Union
import numpy as np
import pandas as pd
from scipy.constants import nano, speed_of_light, pi
from pathlib import Path
import json
from datetime import datetime
import re
import os

time_stamp = datetime.now().strftime("%Y%m%d%H%M")
HOME_PATH = os.environ['HOME']
CALC_STATUS = Path(HOME_PATH)/'projects'/'rrg-aspuru'/'cicsy'/'madness'/'calc_status.csv'

DESCRIPTORS_ABS = ("osc_strength", "electronic_gap", "ext_spectra_com", "ext_peak", "ext_peak_w")
DESCRIPTORS_EM = ("osc_strength", "electronic_gap", "fluo_spectra_com", "fluo_peak_w", "fluo_rate")

TRANSITIONS = {
	"abs_0_1": DESCRIPTORS_ABS,
	"abs_0_2": DESCRIPTORS_ABS,
	"abs_0_3": DESCRIPTORS_ABS,
	"fluo_1_0": DESCRIPTORS_EM
}

def extract_ES_descriptors(es_file):
	es_properties = {}
	with open(es_file, 'r') as f_in:
		for i, line in enumerate(f_in):
			if 0 <= i <5:
				es_properties[f'sr_index_{i+1}'] = float(line.split()[-2])
			if 6 <= i < 11:
				es_properties[f'd_index_{i-5}_x'] = float(line.split()[2])
				es_properties[f'd_index_{i-5}_y'] = float(line.split()[4])
				es_properties[f'd_index_{i-5}_z'] = float(line.split()[6])
				es_properties[f'd_index_{i-5}_sum'] = float(line.split()[9])
			if 12 <= i < 17:
				es_properties[f'hole_rmsd_{i-11}_x'] = float(line.split()[6])
				es_properties[f'hole_rmsd_{i-11}_y'] = float(line.split()[7])
				es_properties[f'hole_rmsd_{i-11}_z'] = float(line.split()[8])
				es_properties[f'hole_rmsd_{i-11}_n'] = float(line.split()[10])
			if 18 <= i < 23:
				es_properties[f'ele_rmsd_{i-17}_x'] = float(line.split()[6])
				es_properties[f'ele_rmsd_{i-17}_y'] = float(line.split()[7])
				es_properties[f'ele_rmsd_{i-17}_z'] = float(line.split()[8])
				es_properties[f'ele_rmsd_{i-17}_n'] = float(line.split()[10])
			if 24 <= i < 29:
				es_properties[f'h_index_{i-23}_x'] = float(line.split()[2])
				es_properties[f'h_index_{i-23}_y'] = float(line.split()[4])
				es_properties[f'h_index_{i-23}_z'] = float(line.split()[6])
				es_properties[f'h_index_{i-23}_ct'] = float(line.split()[8])
				es_properties[f'h_index_{i-23}_i'] = float(line.split()[11])
			if 30 <= i < 35:
				es_properties[f't_index_{i-29}'] = float([x for x in re.split(r'[:\s*]', line) if x][3])
			if 36 <= i < 41:
				es_properties[f'delta_r_{i-35}'] =	float(line.split()[-2])
			if 42 <= i < 47:
				es_properties[f'es_lambda_{i-41}'] =  float(line.split()[-1])
			if (50 <= i < 55) or (58 <= i < 73):
				dipoles = line.split()
				es_properties[f'trans_dipole_{dipoles[0]}_{dipoles[1]}_x'] = float(dipoles[2])
				es_properties[f'trans_dipole_{dipoles[0]}_{dipoles[1]}_y'] = float(dipoles[3])
				es_properties[f'trans_dipole_{dipoles[0]}_{dipoles[1]}_z'] = float(dipoles[4])
				es_properties[f'trans_dipole_{dipoles[0]}_{dipoles[1]}_diff'] = float(dipoles[5])
				es_properties[f'trans_dipole_{dipoles[0]}_{dipoles[1]}_osc'] = float(dipoles[6])
	return es_properties

def extract_descriptors(calc_dir: Path, transitions: Dict[str, tuple], fmos: int = 0) -> Dict[str, float]:
	"""
	Extract descriptors from the processed results of the QChem TDDFT calculations from Cyrille's code.
	Args:
	 calc_dir: Path to the calculation directory.
	 transitions: Dictionary of transitions to analyze, and which parameters to extract for each transition
				  (requires the existence of a "TRANSITION_NAME.npz" file)
	fmos: Number of frontier molecular orbitals to include (0: None, 1: only HOMO and LUMO, 2: also HOMO-1 / LUMO+1, ...)
	Returns:
	Dict[str, float]: Dictionary of descriptor names and values.
	"""
	descriptors: dict = {}

	for transition, desc_keys in transitions.items():
		loaded_data = np.load(calc_dir / f"{transition}.npz")
		for descriptor in desc_keys:
			descriptors[f"{transition}_{descriptor}"] = float(loaded_data[descriptor])
	descriptors["max_gain_wl"], descriptors["max_gain_factor"] = calculate_gain_factor(calc_dir / "specs.csv")
	fmo_energies = orca_parse_orbital_energies(calc_dir / "orca_mo.out") if fmos > 0 else None
	for fmo_idx in range(fmos):
		descriptors[f"homo_-{fmo_idx}"] = fmo_energies[0][-(fmo_idx + 1)]
		descriptors[f"lumo_+{fmo_idx}"] = fmo_energies[1][fmo_idx]

	return descriptors

def orca_parse_orbital_energies(orca_mo_file):
	occ = []
	emp = []
	with open(orca_mo_file, 'r') as f_in:
		for i, line in enumerate(f_in):
			if 5 <= i :
				if '2' in line.split()[1]:
					occ.append(float(line.split()[-2]))
				else:
					emp.append(float(line.split()[-2]))
	return occ, emp

def calculate_gain_factor(spec_file: Path, analysis_range: tuple = (300, 800), abs_threshold: float = 0.05) -> Tuple[float, float]:
	"""
	Takes the calculated spectra (absorption and emission) from Cyrille's code and extracts a gain factor.

	Args:
	spec_file: Path to the specs.csv file
	analysis_range: Tuple (min, max) of the wavelength range to consider for gain factor analysis.
	abs_threshold: Threshold for absorption when determining the maximum gain factor.
	Returns:
	float: Maximum gain wavelength (in nm).
	float: Maximum gain factor (in m^2 s).
	"""
	spectra: pd.DataFrame = pd.read_csv(spec_file, usecols=[1, 2, 3])
	wavelengths: np.array = spectra["x_nm"].to_numpy()
	absorption_norm: np.array = spectra["extinct"].to_numpy() / np.amax(spectra["extinct"].to_numpy())
	emission: np.array = spectra["fluo"].to_numpy()
	
	analysis_idx: np.array = (wavelengths > analysis_range[0]) & (wavelengths < analysis_range[1])
	wavelengths = wavelengths[analysis_idx]
	absorption_norm = absorption_norm[analysis_idx]
	emission = emission[analysis_idx]

	gain_spectrum: np.array = ((wavelengths*nano)**4 * emission) / (8 * pi * speed_of_light * np.trapz(np.flip(emission), np.flip(wavelengths*nano)))

	no_absorption_idx: np.array = absorption_norm < abs_threshold
	max_gain_factor: float = np.amax(gain_spectrum[no_absorption_idx])
	max_gain_wl: float = wavelengths[no_absorption_idx][np.argmax(gain_spectrum[no_absorption_idx])]

	return max_gain_wl, max_gain_factor

def get_smiles(job_file: Path):
	with open(job_file, 'r') as f_in:
		for i, line in enumerate(f_in):
			if 'smi=' in line:
				smi = line[5:-2]
			if 'name=' in line:
				hid = line[6:-2]
	return {'hid': hid, 'smiles': smi}


work_dir = Path.cwd()
print(work_dir)
job_file = work_dir/'madness_job.sh'
#print(job_file, job_file.exists())
es_file = work_dir/'spec'/'ES_results.out'
#print(es_file, es_file.exists())
spec_file = work_dir/'spec'/'specs.csv'
#print(spec_file, spec_file.exists())
spec_dir = work_dir/'spec'
hid = str(work_dir.name)
#print(hid)
if es_file.exists() and spec_file.exists():
	hidsmi = get_smiles(job_file)
	#print(hidsmi.keys(), "rtn0 done")
	rtn1 = extract_ES_descriptors(es_file)
	#print(rtn1.keys(), "rtn1 done")
	rtn2 = extract_descriptors(calc_dir = spec_dir, transitions = TRANSITIONS, fmos= 10)
	#print(rtn2.keys(), 'felix_des')
	des = {**rtn1, **rtn2}
	#print(des.keys(), 'rtn done')
	rtn = { 'hidsmi': hidsmi, 'descriptors': des}
	json_out = work_dir/f'{hid}.json'
	#print(str(json_out))
	with open(json_out, 'w') as f_out:
		json.dump(rtn, f_out)


