import os, sys
import datetime
from Auto_opt_measurements_ocean import Optical_measurements, power_on, power_off


## settings
ml = 1e-3

#initialize the device############
# power_on(['AbsPL', 'TE', 'Evap'])

Opt = Optical_measurements(TE = True, AbsPL = True, Pump = True, Evap=True, logger=True, power_control=False)


# Opt.open_valves() #initialze valve at A and 1

# # self.log('washing selector to pump line with solvent {timestamp}')
# Opt.draw_and_dispense('ACN', 'waste', 4*ml,  draw_velocity = 1000, dispense_velocity = 1000)  
# Opt.log('emptying selector to pump line with the air {timestamp}')
# Opt.draw_and_dispense('air', 'valve', 1*ml,  draw_velocity = 1000, dispense_velocity = 1000)  

# ###turn lamp and LED off#########
Opt.close()

# ###shutdown######################
Opt.__del__()
power_off(['AbsPL', 'TE', 'Evap'])


# Opt.draw_and_dispense('vial', 'waste', 1*ml,  1000, 1000)
# Opt.draw_and_dispense('vial', 'waste',  1*ml, 1000, 1000)
# Opt.vial_wash(0.5*ml, 3)
# Opt.change_solvent(['TE', 'absorption', 'PL'])
# Opt.collector_wash(0.7*ml, 2, 2)
# Opt.collector_wash(0.7*ml, 2, 3)
# Opt.collector_wash(0.7*ml, 2, 4)
# Opt.collector_wash(0.7*ml, 2, 5)
# Opt.collector_wash(0.7*ml, 2, 7)
# Opt.collector_wash(0.7*ml, 2, 8)
# Opt.collector_wash(0.7*ml, 2, 9)
# Opt.collector_wash(0.7*ml, 2, 10)

###test Evap##################

# Opt.do_evaporation([10], evap_temp = 50, duration = 30, after_temp= 35, nitrogen=True, wait=True)
# Opt.do_dissolution([10], solvent_volume=0.15*ml, temparature=35, duration = 300, mixing_speed=8, after_temp=25, wash_line = False, wait = True)

# Opt.open_valves()
# Opt.selector.move_to(10)
# Opt.draw_and_dispense('ACN', 'waste', 5*ml, 1000, 1000)
# Opt.close_valves()

# sample_info ={
#     'name' : 'Blank',
#     'concentration(uM)' : 0,
#     'solvent' : 'ACN'
# }

# path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20210210_BSBCz_toluene'
# if os.path.exists(path) == False:
#     os.mkdir(path)
# fname = '%s/%s_%s_in_%s_1' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])
# fname = '%s/%s_%s_in_%s_with_Evap_after_redislv_1' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])

# # Opt.open_valves()
# # Opt.selector.move_to(10)
# # Opt.draw_and_dispense('ACN', 'valve', 0.2*ml, draw_velocity=500, dispense_velocity=500)
# # Opt.close_valves()

# # Opt.do_redissolution([10])
# Opt.do_measurements(fname, measurements = ['absorption', 'PL', 'TE'], sample_position = None, sample_volume = 0.12*ml, \
#                                                     dilution = 1, solvent = 'toluene', water_ratio = 0, sample_info=sample_info) 
# 
# Opt.measure_blank(fname, measurements = ['absorption', 'PL', 'TE'], sample_position = 2, sample_volume = 0.2*ml, solvent='ACN', sample_info=sample_info)

# Opt.vial_wash(0.3*ml, 4)

####Automeasurement###########
# Opt.auto_measurement(measurements = ['absorption','PL', 'TE'], solvent='ACN', redissolution=False, measure_blank=False)


###manual measurment_abs_PL_TE#############
# sample_info ={
#     'name' : 'BSBCz',
#     'concentration(uM)' : 0,
#     'solvent' : 'toluene'
# }

# path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20210208_Evap_test_toluene'
# if os.path.exists(path) == False:
#     os.mkdir(path)
# fname = '%s/%s_%s_in_%s_filteredx10_nofilter_3' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])
# # fname2 = '%s/%s_%s_in_%s_3_2' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])
# # fname3 = '%s/%s_%s_in_%s_3_3' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])

# Opt.do_measurements(fname, measurements = ['absorption','PL','TE'], sample_position = None, sample_volume = 0.12*ml, \
#                                 solvent = 'toluene', dilution = 1, water_ratio = 0, sample_info=sample_info) 
# Opt.do_measurements(fname2, measurements = ['absorption', 'PL','TE'], sample_position = 3, sample_volume = 0.12*ml, dilution = 1, water_ratio = 0, sample_info=sample_info) 
# Opt.do_measurements(fname3, measurements = ['absorption', 'PL','TE'], sample_position = 4, sample_volume = 0.12*ml, dilution = 1, water_ratio = 0, sample_info=sample_info) 

# Opt.measure_blank(fname, measurements = ['absorption', 'PL','TE'], sample_position = 2, sample_volume = 0.15*ml, sample_info=sample_info) 
# Opt.measure_blank(fname2, measurements = ['absorption', 'PL'], sample_position = 3, sample_volume = 0.2*ml, sample_info=sample_info) 
# Opt.measure_blank(fname3, measurements = ['absorption', 'PL'], sample_position = 3, sample_volume = 0.2*ml, sample_info=sample_info) 
############

# Opt.draw_and_dispense('ACN', 'waste', 5*ml, draw_velocity=1000, dispense_velocity=1000)
# Opt.draw_and_dispense('ACN', 'flow_cell', 2*ml, draw_velocity=500, dispense_velocity=500)
# Opt.draw_and_dispense('ACN', 'flow_cell_PL', 2*ml, draw_velocity=500, dispense_velocity=500)
# # Opt.draw_and_dispense('ACN', 'flow_cell_PL', 0.5*ml, draw_velocity=500, dispense_velocity=500)
# # # Opt.draw_and_dispense('ACN', 'flow_cell_PL', 0.5*ml, draw_velocity=500, dispense_velocity=500)
# Opt.draw_and_dispense('ACN', 'flow_cell_abs', 2*ml, draw_velocity=500, dispense_velocity=500)
# Opt.vial_wash(0.3*ml, 3)
# Opt.collector_wash(0.5*ml, 3, 9)
# Opt.collector_wash(0.5*ml, 3, 10)
# # # Opt.collector_wash(0.3*ml, 1, 8)

# # # # # Opt.selector = Valco_valve('visa://10.22.1.20/ASRL17::INSTR', dev_id = 0, mode = 3, position = 1)
# # # # # Opt.valve = Valco_valve('visa://10.22.1.20/ASRL4::INSTR', dev_id = 0, mode = 1, position = 'A')
# # # # # # Opt.selector.move_to(2)
# # # # # # Opt.draw_and_dispense('ACN', 'valve', 0.2*ml, draw_velocity=500, dispense_velocity=500)
############


#fill cells########
#Opt.pump.draw(volume = 1*ml, valve='valve')
#time.sleep(10)
#Opt.pump.dispense(volume = 1*ml, valve='vial')
#Opt.draw_and_dispense('ACN', 'waste', 1*ml, draw_velocity=500, dispense_velocity=500)
#Opt.draw_and_dispense('ACN', 'waste', 1*ml, draw_velocity=500, dispense_velocity=500)
#Opt.draw_and_dispense('ACN', 'waste', 1*ml, draw_velocity=500, dispense_velocity=500)
#Opt.draw_and_dispense('ACN', 'flow_cell', 0.5*ml, draw_velocity=500, dispense_velocity=500)
#Opt.draw_and_dispense('air', 'valve', 1*ml, draw_velocity=2000, dispense_velocity=2000)
# Opt.draw_and_dispense('air', 'valve', 1*ml, draw_velocity=1000, dispense_velocity=1000)



#fill vial#################
# Opt.selector = Valco_valve('visa://10.22.1.20/ASRL17::INSTR', dev_id = 0, mode = 3, position = 1)
# Opt.valve = Valco_valve('visa://10.22.1.20/ASRL4::INSTR', dev_id = 0, mode = 1, position = 'A')
# Opt.selector.move_to(3)
# Opt.draw_and_dispense('ACN', 'valve', 0.2*ml, draw_velocity=500, dispense_velocity=500)
# # Opt.selector.move_to(10)7
# # Opt.draw_and_dispense('ACN', 'valve', 0.2*ml, draw_velocity=500, dispense_velocity=500)
# Opt.selector.close_device()
# # Opt.valve.close_device()
#time.sleep(10)


#Opt.measure_PL_manual(fname, sample_info=sample_info) 



###shut down########
#
# 
#Opt.close()

###measure line#############
#Opt = Optical_measurements(TE=False, AbsPL=True)

#Opt.mgr.pump.draw_and_dispense('air', 'flow_cell_PL', 1*ml, velocity=1000)
#Opt.mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 1*ml, velocity=1000)
#Opt.measure_line('flow_cell_PL', fill='valve')
#Opt.mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.34*ml, velocity=1000)
#input()
#Opt.mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.39*ml, velocity=1000)

# Opt.pump.draw_and_dispense('valve', 'waste', 1*ml, velocity=1000)
# Opt.pump.draw_and_dispense('ACN', 'waste', 0.1*ml, velocity=1000)
# Opt.measure_line('valve', fill='ACN', volume=0.1*ml, step=0.005*ml)

###########################
