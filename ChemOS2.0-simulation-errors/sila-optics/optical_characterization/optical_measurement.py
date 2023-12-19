import os, sys
import datetime, time
from optical_characterization.Auto_opt_measurements_ocean import Optical_measurements, power_on, power_off, lamp_on, lamp_off, led_on, led_off

"""
This is the script to excecute the measurements for optical characterizations.
 
"""

## settings
ml = 1e-3


"""
Commands to turn on the devices. 
The devices should be turned on before using lamp_on and led_on commands.

power_on : turn on the power of devices.
    parameters : list of devices('AbsPL', 'TE', 'Evap')

lamp_on : turn on the deuterium and halogen lamp (DH-mini). 
led_on : turn on the LED.
    parameters : led power (0-100)
# """
# power_on(['AbsPL','TE'])
# lamp_on()
# led_on(percentage= 30)


"""
Create the instance of optical measurements

parameters :
    TE : should be True if you want to do the PL lifetime measurements.
    AbsPL : should be True if you want to do the absorption and PL measurements.
    Pump : should be True to use the syringe pump (basically need for all the measurements).
    Evap : should be True if you want to do the solvent evaporation and redissolution.
    Logger
    power_control
"""
Opt = Optical_measurements(TE = True, AbsPL = True, Pump = True, Evap=False, logger=True, power_control=False)

"""
Commands for Manual measurement. The measurement parameters are specified in "config_AutoOpt_ocean.py".


Opt.do_redissolution : Evaporate the solvent and redissolute the sample in the solvent for the measurements.
                       To do this, "Evap" should be True (parameters are in  "config_AutoOpt_ocean.py").
    parameters 
        collector_nums : list of collector vial numbers where the sample to be redissolution are (2,3,4,5,7,8,9,10). 

Opt.do_measurements : excecute the measurements
    parameters : 
        fname : filename used to save the results
        measurements : list of the measurements you do ('absorption', 'PL', 'TE').
        sample_position : collection vial number where your sample are (2,3,4,5,7,8,9,10). 
                          The setup will collect sample from the dilution vial if None.
        sample_volume : volume of the sample (e.g. 0.2*ml)
        dilution : put the number larger than 1 if you want to make dilution before the measurements.
                   (e.g. if you set 2, the setup makes 2 times dilution. set 1 if you don't do dilution.)
        solvent : solvent used for the measurements. This is used to calculate the PLQY using its refractive index.
                  currently we have 'ACN', 'toluene' 'Ethanol' in the list.
                  if you want to use new solvent, please add it to the list in "config_AutoOpt_ocean.py".
        water_ratio : to adjest the water contents of the reference solvents. set 0 if the sample solution doesn't contain water.
        fill_ref : send solvent to the flow cells before the measurement if True.
        sample_info : dictionary of sample information for record.

Opt.measure_blank : excecute the blank measurement by using a solvent as blank sample.
    parameters :
        fname : same as above
        measurements : same as above
        sample_position : same as above
        sample_volume : this volume of the solvent is first send to the specified collection vial and start blank measurements.
        solvent : same as above
"""

# sample_info ={'name' : 'Blank', 'concentration(uM)' : 0, 'solvent' : 'ACN'}

# path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20210316_test'
# if os.path.exists(path) == False:
#     os.mkdir(path)
# fname = '%s/%suM_%s_in_%s_3' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])

# Opt.do_measurements(fname, measurements = ['absorption', 'PL', 'TE'], sample_position = None,\
#         sample_volume = 0.12*ml, dilution = 1, solvent = 'ACN', water_ratio = 0, fill_ref = False, sample_info=sample_info) 

# Opt.do_redissolution([2])

# # # Opt.do_measurements(fname, measurements = ['absorption', 'PL', 'TE'], sample_position = None,\
# # #         sample_volume = 0.12*ml, dilution = 1, solvent = 'ACN', water_ratio = 0, fill_ref = False, sample_info=sample_info) 

# Opt.measure_blank(fname, measurements = ['absorption','PL', 'TE'], sample_position = 2, \
#           sample_volume = 0.2*ml, solvent='ACN', fill_ref = False, sample_info=sample_info)


"""
Commands for automatic measurement. The measurement parameters are also specified in "config_AutoOpt_ocean.py".

Opt.auto_measurement : Excecute the automatic measurements. The code waits for a job input from HPLCMS.
    parameters : 
        measurements : list of the measurements you do ('absorption', 'PL', 'TE').
        solvent : solvent used for the measurements. This is used to calculate the PLQY using its refractive index.
                  currently we have 'ACN', 'toluene' 'Ethanol' in the list.
                  if you want to use new solvent, please add it to the list in "config_AutoOpt_ocean.py".
        redissolution : do sample evaporation and redissolution if True. Parameters are in "config_AutoOpt_ocean.py".
        measure_blank : do blank measurements before waiting the job input, if True
# """
# Opt.auto_measurement(measurements = ['absorption','PL', 'TE'], 
#                              solvent='ACN', redissolution=False, measure_blank=False)


"""
Commands for the maintainance.

Opt.vial_wash : wash dilution vial with the solvent. 
    parameters 
        volume : solvent volume to be used for washing.
        repeat : how many times you wash the vial.
        fill : vial is filled with the solvent after the wash, if True.
Opt.collector_wash :  wash collector vials with the solvent. 
    parameters 
        volume : same as above
        repeat : same as above
        collector_nums : list of collector vial numbers to be washed (2,3,4,5,7,8,9,10)
        fill : vials are filled with the solvent after the wash, if True.
Opt.cell_wash :   wash measurement cells with the solvent.   
    parameters    
        cells : list of flow cells to be washed ('flow_cell', 'flow_cell_abs', 'flow_cell_PL').
        volume : same as above
        repeat : same as above
Opt.vial_empty :  send solution in the dilution vial to waste.
Opt.collector_empty : send solutions in the collector vials to waste.
    parameters 
        collector_nums : list of collector vial numbers to be washed (2,3,4,5,7,8,9,10)
# """
Opt.cell_wash(['flow_cell', 'flow_cell_abs', 'flow_cell_PL'], volume = 1*ml, repeat = 2)
Opt.vial_wash(0.3*ml, repeat=2, fill = False)
Opt.collector_wash(0.5*ml, repeat=1, collector_nums = [2,3,4,5,7,8,9,10], fill = False)
 
# Opt.vial_empty()
# Opt.collector_empty([2,3])

"""
Commands to turn off the devices.
Instance of the optical measurements should be deleted before using the power_off command.

Opt.close : turn off lamps and led
Opt.__del__ : delete the instance of optical measurements
power_off : turn off the power of devices.
    parameters : list of devices('AbsPL', 'TE', 'Evap')
"""
###turn lamp and LED off#########
Opt.close()

###shutdown######################
Opt.__del__()
power_off(['AbsPL', 'TE', 'Evap'])

"""
Commands used to modify the setup.

Opt.measure_line :
    parameters
        line : line to be measured ('flow_cell', 'flow_cell_abs', 'flow_cell_PL', 'vial', 'valve')
        direction : 'dispense' or 'draw'
        volume : volume of the solvent used to measure line volume (default 1*ml)
        step : step of the solvent to draw from or dispense to the line (default 0.01*ml)
"""

# Opt.measure_line('flow_cell_abs', direction = 'dispense', volume = 0.5*ml, step=0.01*ml)











###test Evap##################

# Opt.do_evaporation([10], evap_temp = 50, duration = 30, after_temp= 35, nitrogen=True, wait=True)
# Opt.do_dissolution([10], solvent_volume=0.15*ml, temparature=35, duration = 300, mixing_speed=8, after_temp=25, wash_line = False, wait = True)
# Opt.collector_wash(0.3*ml, 4, [2])
# Opt.open_valves()
# Opt.selector.move_to(10)
# Opt.draw_and_dispense('ACN', 'waste', 5*ml, 1000, 1000)
# Opt.close_valves()

###########################
