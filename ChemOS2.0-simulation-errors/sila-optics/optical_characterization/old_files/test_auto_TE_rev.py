from pylab.manager import Manager, Logger
import pylab.instruments as ins
import numpy as np
import matplotlib.pyplot as plt
import time
from transient_emission import Transient_Emission

## settings
ml = 1e-3
loop_to_pump  = 0.04*ml # (under)
vial_to_pump = 0.15*ml
to_cell_under = 0.46*ml #0.43
to_cell_over  = 0.51*ml #0.57


# trans_config = 75, 10e-3 # about 3 (add OD2 to the white light lamp)
# pl_config    = 8, 0.2 # about *2
# trans_config = 2000, 10e-3 # about *3 (add OD2 to the white light lamp)
# pl_config    = 150, 0.2 # about *2

log = Logger(stdout=True, logfile=None, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')
mgr = Manager([

    {
        "name": "pump",
        "instrument": "PumpPSD8",
        "settings": {
            "visa": "ASRL2::INSTR",
            "addr": 0,
            "syringe_volume": 1E-3,
            "init_valve": 8,
            "ports": {
                "vial": 1,
                "ACN" : 2,
                "air": 3,
                "flow_cell": 5,
                "waste": 7 ,
                "valve": 8,
            }
        }
    },

])
mgr.pump.set_velocity(2000)

## measurement functions
delay = lambda t=0.1: time.sleep(t)

def measure_line(line, fill=None):
    if fill:
        mgr.pump.draw_full(fill)
        action = mgr.pump.dispense
    else:
        action = mgr.pump.draw

    c = ''
    volume = 0
    while c != 'c':
        action(volume=0.01*ml, valve=line)
        volume += 0.01
        c = input(f'{volume:.3f}')

    mgr.pump.dispense_all('waste')

def save_data(fname, data):
    np.savetxt(fname, data, delimiter=',', header='wl (nm), abs (%), pl (and abs/pl pair repeats)')

def vial_wash(volume):
    for _ in range(3):
        mgr.pump.draw_and_dispense('valve', 'vial', volume, velocity=1000)
        mgr.pump.draw_and_dispense('vial', 'waste', volume*1.2, velocity=2000)

## start pl and abs measurement with pump
def measure_TE(fname, dilution, TE_setting):
    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    mgr.valve._move_to('A')
    _ = input('Please load your material')

    # switch position
    log('Switch valve to inject position')
    mgr.valve._move_to('B')

    # draw sample and extra THF for dilutions
    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=200)

    total_volume = dilution*0.1*ml
    print(total_volume)
    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=200)
    mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity =1000) # mix 3 times


    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    flow_cell_volume = to_cell_over - to_cell_under + 2*buffer_volume
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')
    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    mgr.pump.dispense_all(valve='flow_cell')
    mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_under-buffer_volume, velocity=1000)

    log('Measure transient emission {timestamp}')
    TE.measure_TE(**TE_setting)

    result = TE.measure_TE(save_filename = fname, name = 'di-nBu-PTPTP', concentration = "500uM", solvent = 'THF', fit_order = 'auto', fit_weight = 0.7, do_plot = True)
    
    # send thf to flow cell and measure background
    log('Wash flow cell')
    for _ in range(4):
        mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_over, velocity=500)

    #clean the flow cell
    log('Wash vial')
    for _ in range(3):
        mgr.pump.draw_and_dispense('valve', 'vial', 1.5*total_volume, velocity=500)
        mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)



def measure_TE_from_vial(fname, TE_setting):
    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    # mgr.valve._move_to('A')
    # _ = input('Please load your material')

    # # switch position
    # log('Switch valve to inject position')
    # mgr.valve._move_to('B')

    # draw sample and extra THF for dilutions
    log('Discard vial to pump volumes')
    mgr.pump.draw_and_dispense('vial', 'waste', vial_to_pump, velocity=200)

    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    flow_cell_volume = to_cell_over - to_cell_under + 2*buffer_volume
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_and_dispense('vial', 'flow_cell', flow_cell_volume)
    mgr.pump.draw_and_dispense('ACN', 'flow_cell', to_cell_under-buffer_volume, velocity=1000)

    log('Measure transient emission {timestamp}')
    TE.measure_TE(**TE_setting)

    TE.measure_TE(save_filename = fname, **TE_setting)
    
    # send thf to flow cell and measure background
    log('Wash flow cell')
    for _ in range(4):
        mgr.pump.draw_and_dispense('ACN', 'flow_cell', to_cell_over, velocity=500)




if __name__ == '__main__':
    

    TE_setting ={
        'name' : 'di-nBu-PTPTP',
        'concentration' : '500uM',
        'solvent' : 'THF',
        'fit_order' : 'auto', 
        'fit_weight' : 0.7, 
        'do_plot' : True
    }

    

    #measure_line('vial', fill=None)
    #measure_line('vial', fill='ACN')
    #mgr.pump.draw_and_dispense('air', 'flow_cell', 3*ml)
    #mgr.pump.draw_and_dispense('air', 'vial', 3*ml)
    measure_TE_from_vial(None, TE_setting)
    #path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200306_test_autoPL/'
    #sample_name = '20uL-di-nBu-PTPTP_in_ACN_dilute1' 
    #measure_TE(path + sample_name, 1, TE_setting)
    # measure_abs_pl('dsb1.csv', dilution=10)
    # measure_abs_pl('ptptp.csv', dilution=5)
    #mgr.pump.draw_and_dispense('air', 'valve', 5*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('air', 'waste', 5*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('valve', 'vial', 0.2*ml, velocity = 1000)
    #mgr.pump.draw(0.5*ml, 'valve')
    #mgr.pump.dispense(0.5*ml, 'vial')
    #mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell', 0.4*ml, velocity = 500)
    #mgr.pump.draw_and_dispense('valve', 'waste', 0.3*ml, velocity = 500)
    #mgr.pump.draw(0.1*ml, 'ACN')
    #mgr.pump.draw('ACN', 0.1*ml, velocity = 500)
    #mgr.pump.dispense(0.04*ml, 'valve')
    
    #vial_wash(1*ml)


    #mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    #measure_line('flow_cell', fill='ACN')
    #measure_line('vial', fill='ACN')
   # mgr.pump.draw_and_dispense('ACN', 'vial', 0.05*ml)
    #mgr.pump.draw_and_dispense('vial', 'waste', 1*ml)
    #mgr.pump.draw_and_dispense('valve', "waste", 0.2*ml, velocity = 500)

    # # extra wash
    # mgr.valve._move_to('B')
    # for _ in range(3):
    #     mgr.pump.draw_and_dispense('valve', 'vial', 1.5*ml, velocity=500)
    #     mgr.pump.draw_and_dispense('vial', 'waste', 2.0*ml, velocity=2000)
    # for _ in range(3):
    #     mgr.pump.draw_and_dispense('valve', 'flow_cell', 1.0*ml, velocity=500)
    # mgr.pump.draw_and_dispense('valve', 'vial', 1*ml, velocity=500)

    # pass






