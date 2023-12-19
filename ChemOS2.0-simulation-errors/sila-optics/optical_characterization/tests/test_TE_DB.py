from transient_emission import Transient_Emission
import matplotlib.pyplot as plt
from pprint import pprint as pp

TE = Transient_Emission(device = True, DB = False)

###measurement##########
"""
experimental parameters are in ".config_transient_emission.dat" 
"""
TE.detector_on()

sample_info ={
    'name' : 'di-nBu-PTPTP',
    'concentration' : '500uM',
    'solvent' : 'THF'
}

result = TE.measure_TE(save_filename = None, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **sample_info)

#syn_id = TE.db.get('synthesis')['synthesis_id'][0]
#TE.db.add_TE_results(result, synthesis_id=syn_id)


####enumrator#######################
#data = TE.load_pkl('C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/for_automation/sample_measured/2uL-500uM-di-nBu-PTPTP_ACN100_char4_di-nBu-PTPTP_PL_1/2uL-500uM-di-nBu-PTPTP_ACN100_char4_di-nBu-PTPTP_PL.pkl')
#data = TE.load_pkl('C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/for_automation/sample_measured/2uL-500uM-diphenylanthracene_ACN100_char_diphenylanthracene_PL/2uL-500uM-diphenylanthracene_ACN100_char_diphenylanthracene_PL.pkl')
##data = TE.load_pkl('C:/Users/smily/Dropbox/PythonScript/kazu/data/20200227_optical_table_test/2020027_optical_table_test_IRF.pkl')
#
#fname = 'D:/test'
#result = TE.measure_TE_emu(data['raw'], 2E7, save_filename = None, name = 'di-nBu-PTPTP', concentration = "500uM", solvent = 'THF', fit_order = 'auto', fit_weight = 0.7, do_plot = True)

###DB###############################

# experiments = TE.db.get_local_experiments()
# print(experiments)
# id = experiments['experiment_id'][1]

# print(TE.db.get_experiment(id)['metadata'][0])
# TE.db.plot_result(id)
# df = TE.db.get('xy_data', filters =[TE.db.client.models.xy_data.experiment_id == id])
# print(TE.db.get_raw_data(id)[0])

# TE.db._apply_config()
