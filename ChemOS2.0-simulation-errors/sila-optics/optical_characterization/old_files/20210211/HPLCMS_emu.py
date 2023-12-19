import pickle as pkl

folder = 'C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/sample_to_measure'

fname = '%s/blank.pkl' %folder
fname = '%s/shutdown.pkl' %folder

data = {}

with open (fname, 'wb') as f:
    pkl.dump(data, f)