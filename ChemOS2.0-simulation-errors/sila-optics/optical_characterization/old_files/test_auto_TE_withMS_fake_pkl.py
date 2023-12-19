import pickle

filename = 'C:/Users/Thermo/Dropbox/PythonScript/kazu/data/for_automation/sample_to_measure/%s_%s.pkl' %('a','b')

characterization_params = {
                            'injection_name' : 'fake',
                            'target_name' : 'di-nBu-PTPTP',
                            'average_absorbance_peak' : 30,
                            'average_absorbance_375' : 20,
                            'loop_volume' : 100
                            }


with open(filename, 'wb') as f:
    pickle.dump(characterization_params, f)
