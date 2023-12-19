from pylab.manager import Manager, Logger
import pylab.instruments as ins
import numpy as np
import matplotlib.pyplot as plt
import time
from transient_emission import Transient_Emission

## settings
ml = 1e-3
loop_to_pump  = 0.04*ml # (under)
to_cell_under = 0.46*ml #0.43
to_cell_over  = 0.51*ml #0.57

to_cell_under_abs  = 0.46*ml 
to_cell_over_abs  = 0.63*ml

# TE_setting = {'min_rate' : 1E-3, 
#               'max_rate' : 3E-2,
#               'accumulation' : 10000, 
#                'filename' : None, 
#                'fit_param' : 0.5 ,
#                'filter_init_position' : 12, 
#                'initial_freq' :  1e6,
#                'do_plot' : False
#               }


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
                "flow_cell_abs": 4,
                "flow_cell": 5,
                "waste": 7 ,
                "valve": 8,
            }
        }
    },

    # {
    #     "name": "spectrometer",
    #     "instrument": "ThorlabsCCS",
    #     "settings": {"visa": "USB0::0x1313::0x8089::M00562925::RAW"}
    # },

    # {
    #     "name": "lamp",
    #     "instrument": "DH_MINI",
    #     "settings": {"visa": 'ASRL8::INSTR',
    #                  "verbose" : True
    #     }
    # },    


    # {
    #     "name": "valve",
    #     "instrument": "Valco_valve",
    #     "settings": {
    #         "visa": "ASRL4::INSTR",
    #         "dev_id" :0, 
    #         "mode" : 1, 
    #         "position" : 'B' 
    #     }
    # }

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
def measure_TE(fname, dilution, **kwargs):

    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    #mgr.valve.move_to('A')
    _ = input('Please load your material')

    # switch position
    log('Switch valve to inject position')
    #mgr.valve.move_to('B')

    # draw sample and extra THF for dilutions

    air_gap_volume = 0.05*ml

    mgr.pump.draw_and_dispense('air', 'flow_cell', air_gap_volume, velocity=200)

    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=200)

    total_volume = dilution*0.1*ml

    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)
    mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity =1000) # mix 3 times


    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    flow_cell_volume = to_cell_over - to_cell_under + 2*buffer_volume
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')

    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    # mgr.pump.dispense(volume=1*ml - flow_cell_volume - air_gap_volume, valve='waste')
    mgr.pump.dispense_all(valve='flow_cell')

    mgr.pump.draw_and_dispense('air', 'flow_cell', air_gap_volume, velocity=200)
    mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_under-buffer_volume - air_gap_volume, velocity=300)

    log('Measure transient emission {timestamp}')
    #TE.measure_TE(**kwargs)

    result = TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)

    # send thf to flow cell and measure background
    log('Wash flow cell')
    for _ in range(4):
        mgr.pump.draw_and_dispense('ACN', 'flow_cell', to_cell_over, velocity=500)
    for _ in range(4):
        mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_over, velocity=1000)

    #clean the flow cell
    log('Wash vial')
    for _ in range(3):
        mgr.pump.draw_and_dispense('valve', 'vial', 3*total_volume, velocity=1000)
        mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    mgr.pump.draw_and_dispense('vial', 'waste', 0.5*ml, velocity=1000)


    #try:
    #_ = TE.measure_TE(save_filename = None, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)
    #except:
     #   pass

if __name__ == '__main__':
    
 #   vial_wash(1*ml)
#     mgr.pump.draw_and_dispense('vial', 'waste', 2*ml, velocity = 1000)
# #    mgr.pump.draw_and_dispense('valve', 'flow_cell', 2*ml, velocity = 500)

    ###vial wash####
    # TE = Transient_Emission(device = True, DB = False)
    # for i in range(3):
    #     mgr.pump.draw_and_dispense('ACN', 'vial', 0.5*ml, velocity=1000)
    #     mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    # for i in range(3):
    #     mgr.pump.draw_and_dispense('valve', 'vial', 0.5*ml, velocity=1000)
    #     mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    # mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    # #####
    

    ###flow cell wash###
    # TE = Transient_Emission(device = True, DB = False)
    # for i in range(5):
    #     mgr.pump.draw_and_dispense('ACN', 'flow_cell', 0.5*ml, velocity=500)
    # for i in range(5):
    #     mgr.pump.draw_and_dispense('valve', 'flow_cell', 0.5*ml, velocity=1000)
    #####


    ###measure###
    conc = '0uM'

    sample_info ={
        'name' : 'diphenylanthracene',
        # 'name' : 'di-nBu-PTPTP',
        'concentration' : conc,
        'solvent' : 'ACN'
    }

    import os
    path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200727_test_autoPL/'
    if os.path.exists(path) == False:
        os.mkdir(path)
    sample_name = '%s-%s_in_ACN' %(conc, sample_info['name']) 
    measure_TE(path + sample_name, dilution = 1, **sample_info)
    ######


    #mgr.pump.draw_and_dispense('air', 'valve', 5*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('air', 'waste', 5*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('valve', 'waste', 1*ml, velocity = 1000)
    
    #mgr.pump.draw_and_dispense('valve', 'vial', 0.2*ml, velocity = 1000)
    #mgr.pump.draw(0.5*ml, 'valve')
    #mgr.pump.dispense(0.5*ml, 'vial')
    #mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell', 0.4*ml, velocity = 500)
    #mgr.pump.draw_and_dispense('valve', 'flow_cell', 2*ml, velocity = 500)
    #mgr.pump.draw_and_dispense('valve', 'waste', 0.3*ml, velocity = 500)
    #mgr.pump.draw(0.1*ml, 'ACN')
    #mgr.pump.draw('ACN', 0.1*ml, velocity = 500)
    #mgr.pump.dispense(0.04*ml, 'valve')
    
    #vial_wash(1*ml)


    #mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    #measure_line('flow_cell_abs', fill='ACN')
    #mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 2*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', 0.63*ml, velocity = 1000)

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






