from pylab.manager import Manager, Logger
import pylab.instruments as ins
import numpy as np
import matplotlib.pyplot as plt
import time
from transient_emission import Transient_Emission
from absorption_and_PL import Abs_PL

## settings
ml = 1e-3
loop_to_pump  = 0.04*ml # (under)
to_cell_under = 0.46*ml #0.43
to_cell_over  = 0.51*ml #0.57

to_cell_under_abs  = 0.35*ml 
to_cell_over_abs  = 0.50*ml
# to_cell_under_abs  = 0.46*ml 
# to_cell_over_abs  = 0.63*ml


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
mgr.pump.set_velocity(1000)

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

    mgr.valve.move_to('A')
    _ = input('Please load your material')

    # switch position
    log('Switch valve to inject position')
    mgr.valve.move_to('B')

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
    #TE.measure_TE(**kwargs)

    result = TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)

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


## start pl and abs measurement with pump
def measure_TE_air(fname, dilution, **kwargs):
    #the flow cell shoudl be filled with solvent at the initial state

    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    mgr.valve.move_to('A')
    _ = input('Please load your material')

    # switch position
    log('Switch valve to inject position')
    mgr.valve.move_to('B')

    log('add air gap to the flow cell')
    mgr.pump.draw_and_dispense('air', 'flow_cell', 0.01*ml, velocity = 200)

    # draw sample and extra THF for dilutions
    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=200)

    total_volume = dilution*0.1*ml
    print(total_volume)
    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=200)
    mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity =1000) # mix 3 times

    #mgr.pump.draw_and_dispense('air', 'vial', 1*ml, velocity =1000) #to make sure all the sample are go to the vial 


    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    flow_cell_volume = to_cell_over - to_cell_under + 2*buffer_volume
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')
    #mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    mgr.pump.dispense(volume=1*ml - total_volume, valve='waste')  #discard air (total sample volume should be less than 1 ml)

    if total_volume > to_cell_over + buffer_volume:
        mgr.pump.dispense(volume = to_cell_over + buffer_volume, valve = 'flow_cell', velocity=1000)
    else:
        mgr.pump.dispense_all(valve='flow_cell')
        #mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_under-buffer_volume, velocity=1000)
        mgr.pump.draw_and_dispense('air', 'flow_cell', to_cell_under-buffer_volume, velocity=1000)

    log('Measure transient emission {timestamp}')
    #TE.measure_TE(**kwargs)

    result = TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)

    log('collect sample to the vial')
    mgr.pump.draw_and_dispense()

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


def measure_abs_PL_manual(fname, dilution, cmd):
    #currently cell volume is 140uL so dilution should be more than 1.5

    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    flow_cell_volume = to_cell_over_abs - to_cell_under_abs + 2*buffer_volume
    air_gap_volume = 0.05*ml

    total_volume = dilution*0.1*ml

    if cmd == 'send_solvent':

        #make sure to fill the flow cell with solvent
        log('Send solvent to the flow cell')
        mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', flow_cell_volume, velocity=200)  ####
        mgr.pump.draw_and_dispense('air', 'flow_cell_abs', air_gap_volume, velocity=200)  ####
        
        # draw sample and extra THF for dilutions
        log('Discard valve to pump volumes')
        mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=200)

        log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
        mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=200)
        mgr.pump.draw_and_dispense('vial', 'vial', 2*ml, velocity =1000) # mix 3 times


        log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
        mgr.pump.draw_full(valve='vial')

        #make air gap between solvent and the sample
        #mgr.pump.dispense(volume = buffer_volume, valve= 'flow_cell')

        mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
        mgr.pump.dispense_all(valve='flow_cell_abs')

        #send solvent to the flow cell
        mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_under - flow_cell_volume - buffer_volume - air_gap_volume, velocity = 200)

    elif cmd == 'send_sample':
    
        mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', flow_cell_volume + air_gap_volume, velocity=200)
        #mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_under_abs-buffer_volume, velocity=1000)

    elif cmd == 'wash':

        # send thf to flow cell and measure background
        log('Wash flow cell')
        for _ in range(4):
            mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_over, velocity=500)

        #clean the flow cell
        log('Wash vial')
        for _ in range(3):
            mgr.pump.draw_and_dispense('valve', 'vial', 1.5*total_volume, velocity=500)
            mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
        mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)


def measure_abs_PL_manual_injection(fname, dark_average = 100, abs_average = 1000, abs_exposure = 0.005, PL_average = 300, PL_exposure = 0.1, filter_size = 20, uv_average = 20):
    #currently cell volume is 140uL so dilution should be more than 2

    results = {}

    #initialize the measurement device
    AP = Abs_PL(led_power=95)

    mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', 0.2*ml, velocity=200)
    
    time.sleep(10)
    input('measure ref')
    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    # log('Measure uv power {timestamp}')
    # uv_ref = AP.measure_uv_power(counts = uv_average, led_power = 95)
    # #time.sleep(5)
    # log('Measure background spectrum {timestamp}')

    # #AP.measure_dark_spectrum(PL_average, PL_exposure, filter_size=filter_size,  do_plot=False)
    # #PL_ref = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    # AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    # abs_ref = AP.measure_transmission_spectrum(abs_average, abs_exposure, filter_size=filter_size, dark_correction=True, do_plot=False)

    mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', 0.2*ml, velocity=200)
    input('measure sample')

    # AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)
    # # time.sleep(10)
    # #input('please load the sample solution')

    # log('Measure absorbance {timestamp}')
    # #PL_ref = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = False, spectral_correction = True, do_plot = False)
    # AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    # abs_sample, transmission, absorption = AP.measure_absorption_spectrum(abs_average, abs_exposure, abs_ref, filter_size=filter_size, dark_correction=True, do_plot=True)
    # #abs_sample = AP.measure_transmission_spectrum(abs_average, abs_exposure, filter_size=filter_size, dark_correction=True, do_plot=False)
    # # plt.plot(abs_sample[0], abs_ref[1])
    # # plt.plot(abs_sample[0], abs_sample[1])
    # #lt.plot(abs_sample[0], abs_sample[1]/abs_ref[1])
    # #plt.xlim(200, 800)
    # #plt.show()

    # log('Measure PL {timestamp}')
    # AP.measure_dark_spectrum(5, PL_exposure, filter_size=filter_size,  do_plot=False)
    # PL_sample = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    # #PL_sample = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = False, spectral_correction = True, do_plot = False)
    # #plt.plot(PL_sample[0], PL_sample[1]-PL_ref[1])
    # plt.plot(PL_sample[0], PL_sample[1])
    # plt.xlim(300, 800)
    # plt.show()

    # log('Measure uv power {timestamp}')
    #uv_sample = AP.measure_uv_power(counts = uv_average, led_power = 95)


    # results = {
    #     'uv' : {
    #         'reference' : uv_ref,
    #         'sample' : uv_sample,
    #         'absorbed_power(W)' : uv_ref - uv_sample 
    #     },
    #     'absorption' : {
    #         'reference' : abs_ref,
    #         'sample' : abs_sample,
    #         'transimittance' : transmission,
    #         'absorbance' : absorption
    #     },
    #     'PL' : PL_sample
    # }

    if fname:
        AP.save_pkl(fname + '.pkl', results)
        np.savetxt(fname + '_abs.csv', np.transpose(np.asarray([abs_ref[0], abs_ref[1], abs_sample[1], transmission[1], absorption[1]])))
        np.savetxt(fname + '_PL.csv', np.transpose(np.asarray([PL_sample[0], PL_sample[1]])))
        # np.savetxt(fname + '_PL.csv', np.transpose(np.asarray(PL_ref[0], PL_ref[1], PL_sample[1], PL_sample[1]-PL_ref[1])))
        np.savetxt(fname + '_uv.csv', np.transpose(np.asarray([uv_ref, uv_sample, uv_ref-uv_sample])))

    return results


## start pl and abs measurement with pump 
def measure_abs_PL(fname, dilution, dark_average = 50, abs_average = 1000, abs_exposure = 0.005, PL_average = 300, PL_exposure = 0.01, filter_size = 20, uv_average = 20):
    #currently cell volume is 140uL so dilution should be more than 2

    results = {}
    
    # calculation to send to flow_cell
    buffer_volume = 0.02*ml
    air_gap_volume = 0.05*ml
    flow_cell_volume = to_cell_over_abs - to_cell_under_abs + 2*buffer_volume

    total_volume = dilution*0.1*ml

    #initialize the measurement device
    AP = Abs_PL(led_power=95)

    #make sure the solvent is in the flowcell
    log('Send solvent to the flow cell')
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', flow_cell_volume, velocity=50)  ####
    #mgr.pump.draw_and_dispense('air', 'flow_cell_abs', air_gap_volume, velocity=200)  ####

    mgr.pump.draw(volume= flow_cell_volume, valve='ACN')
    mgr.pump.set_velocity(200)
    mgr.pump.dispense_all(valve='flow_cell_abs')
    mgr.pump.set_velocity(1000)
    # draw sample and extra solvent for dilutions
    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=200)

    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)

    if dilution > 1:
        for i in range(2):
            mgr.pump.draw_and_dispense('vial', 'vial', total_volume *2, velocity =1000) # mix 3 times
            # mgr.pump.draw_and_dispense('vial', 'vial', 1*ml, velocity =1000) # mix 3 times


    # send sample to the flow cell tubing
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')

    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    mgr.pump.set_velocity(50)
    mgr.pump.dispense_all(valve='flow_cell_abs')
    mgr.pump.set_velocity(1000)

    time.sleep(5)  # need sometime until the solution reaches equiliblium
    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    
    #measure backgrounds
    log('Measure uv power {timestamp}')
    #uv
    uv_ref = AP.measure_uv_power(counts = uv_average, led_power = 95)
    #time.sleep(1)
    log('Measure background spectrum {timestamp}')
    #AP.measure_dark_spectrum(PL_average, PL_exposure, filter_size=filter_size,  do_plot=False)
    #PL_ref = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    #AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    abs_ref = AP.measure_transmission_spectrum(abs_average, abs_exposure, filter_size=filter_size, dark_correction=True, do_plot=False)
    #input('pause')


    #send sample to the flowcell
    #mgr.pump.draw(volume = to_cell_under_abs - buffer_volume, valve = 'valve')
    mgr.pump.draw(volume = to_cell_under_abs - buffer_volume, valve = 'air')
    mgr.pump.set_velocity(50)
    mgr.pump.dispense_all(valve = 'flow_cell_abs')

    # mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_under_abs - buffer_volume, velocity=50)
    #mgr.pump.draw_and_dispense('air', 'flow_cell_abs', air_gap_volume, velocity=200)
    #mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_under-buffer_volume - air_gap_volume, velocity=200)

    time.sleep(15) # need sometime until the solution reaches equiliblium
    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    #measurement
    log('Measure absorbance {timestamp}')
    AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    abs_sample, transmission, absorption = AP.measure_absorption_spectrum(abs_average, abs_exposure, abs_ref, filter_size=filter_size, dark_correction=True, do_plot=True)

    log('Measure PL {timestamp}')
    AP.measure_dark_spectrum(5, PL_exposure, filter_size=filter_size,  do_plot=False)
    PL_sample = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    plt.plot(PL_sample[0], PL_sample[1])
    plt.xlim(300, 800)
    plt.show()

    log('Measure uv power {timestamp}')
    uv_sample = AP.measure_uv_power(counts = uv_average, led_power = 95)

    input('c')

    log('Wash flow cell')
    mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 1*ml, velocity=1000)
    for _ in range(4):
        mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_over, velocity=500)

    #clean the flow cell
    log('Wash vial')
    for _ in range(3):
        mgr.pump.draw_and_dispense('valve', 'vial', 1.5*total_volume, velocity=1000)
        mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    results = {
        'uv' : {
            'reference' : uv_ref,
            'sample' : uv_sample,
            'absorbed_power(W)' : uv_ref - uv_sample 
        },
        'absorption' : {
            'reference' : abs_ref,
            'sample' : abs_sample,
            'transimittance' : transmission,
            'absorbance' : absorption
        },
        'PL' :  PL_sample
        # 'PL' : {
            # 'reference' : PL_ref,
            # 'sample' : PL_sample,
            # 'PL' : np.asarray(PL_sample[0], PL_sample[1] - PL_ref[1])
        # }
    }

    if fname:
        AP.save_pkl(fname + '.pkl', results)
        np.savetxt(fname + '_abs.csv', np.transpose(np.asarray([abs_ref[0], abs_ref[1], abs_sample[1], transmission[1], absorption[1]])))
        np.savetxt(fname + '_PL.csv', np.transpose(np.asarray([PL_sample[0], PL_sample[1]])))
        # np.savetxt(fname + '_PL.csv', np.transpose(np.asarray(PL_ref[0], PL_ref[1], PL_sample[1], PL_sample[1]-PL_ref[1])))
        np.savetxt(fname + '_uv.csv', np.transpose(np.asarray([uv_ref, uv_sample, uv_ref-uv_sample])))

    return results


def measure_abs_PL2(fname, dilution, dark_average = 50, abs_average = 1000, abs_exposure = 0.005, PL_average = 300, PL_exposure = 0.01, led_power = 10, filter_size = 20, uv_average = 20):
    #currently cell volume is 140uL so dilution should be more than 2

    results = {}
    
    # calculation to send to flow_cell
    buffer_volume = 0.025*ml
    air_gap_volume = 0.05*ml
    flow_cell_volume = to_cell_over_abs - to_cell_under_abs + 2*buffer_volume

    total_volume = dilution*0.1*ml

    #initialize the measurement device
    AP = Abs_PL(led_power=led_power)

    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=100)

    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)

    if dilution > 1:
        for i in range(2):
            mgr.pump.draw_and_dispense('vial', 'vial', total_volume *2, velocity =1000) # mix 3 times
            # mgr.pump.draw_and_dispense('vial', 'vial', 1*ml, velocity =1000) # mix 3 times

    #measure backgrounds
    time.sleep(5)  # need sometime until the solution reaches equiliblium
    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    
    log('Measure uv power {timestamp}')
    #uv
    uv_ref = AP.measure_uv_power(counts = uv_average, led_power = 5)
    #time.sleep(1)
    log('Measure background spectrum {timestamp}')
    AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    #AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = 95, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
    abs_ref = AP.measure_transmission_spectrum(abs_average, abs_exposure, filter_size=filter_size, dark_correction=True, do_plot=False)
    #input('pause')

    #fill flow cell with air
    mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 2*ml, velocity=2000)

    # send sample to the flow cell tubing
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')

    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    mgr.pump.set_velocity(1000)
    mgr.pump.dispense_all(valve='flow_cell_abs')
    mgr.pump.set_velocity(1000)

    #send sample to the flowcell
    mgr.pump.draw(volume = to_cell_under_abs - buffer_volume, valve = 'air')
    mgr.pump.set_velocity(1000)
    mgr.pump.dispense_all(valve = 'flow_cell_abs')

    # need sometime until the solution reaches equiliblium
    for i in range(5):
        print('waiting %ss' %(50 - i*10))
        time.sleep(10) 

    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    #measurement
    log('Measure absorbance {timestamp}')
    AP.measure_dark_spectrum(dark_average, abs_exposure, filter_size=filter_size,  do_plot=False)
    abs_sample, transmission, absorption = AP.measure_absorption_spectrum(abs_average, abs_exposure, abs_ref, filter_size=filter_size, dark_correction=True, do_plot=False)


    log('Measure uv power {timestamp}')
    uv_sample_before = AP.measure_uv_power(counts = uv_average, led_power = led_power)

    log('Measure PL {timestamp}')
    AP.measure_dark_spectrum(PL_average, PL_exposure, filter_size=filter_size,  do_plot=False) #dark for PL
    PL_raw = AP.measure_PL_spectrum(PL_average, PL_exposure,  led_power = led_power, filter_size = filter_size, dark_correction = True, spectral_correction = True, do_plot = False)

    log('Measure uv power {timestamp}')
    uv_sample_after = AP.measure_uv_power(counts = uv_average, led_power = led_power)


    #calculation
    uv_abs = uv_ref - (uv_sample_before + uv_sample_after)/2
    PL_energy = np.asarray([PL_raw[0], PL_raw[1]/PL_exposure])
    PL_photons = np.asarray([PL_raw[0], PL_energy[1]*PL_energy[0]*1E-9])
    relative_QE = np.sum(PL_photons[1])/uv_abs
    uv_change = (uv_sample_after-uv_sample_before)/uv_sample_before

    print('uv_change(%) : {:.3f}'.format(uv_change*100))
    print('relative_QE : {:.3f}'.format(relative_QE))

    results = {
        'uv' : {
            'reference' : uv_ref,
            'before' : uv_sample_before,
            'after' : uv_sample_after,
            'change' : uv_change,
            'absorbed_power(W)' : uv_abs 
        },
        'absorption' : {
            'reference' : abs_ref,
            'sample' : abs_sample,
            'transimittance' : transmission,
            'absorbance' : absorption
        },
        'PL' :  {
            'energy' : PL_energy,
            'photons' : PL_photons,
            'relative_QE' : relative_QE
        } 
    }
    
    if fname:
        AP.save_pkl(fname + '.pkl', results)
        np.savetxt(fname + '_abs.csv', np.transpose(np.asarray([abs_ref[0], abs_ref[1], abs_sample[1], transmission[1], absorption[1]])))
        np.savetxt(fname + '_PL.csv', np.transpose(np.asarray([PL_energy[0], PL_energy[1], PL_photons[1]])))
        np.savetxt(fname + '_uv.csv', np.transpose(np.asarray([uv_ref, uv_sample_before, uv_sample_after, uv_change, uv_abs])))
        #absoption spectrum
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
        ax1.plot(absorption[0], absorption[1])
        plt.xlim(300,800)
        plt.savefig('%s_absorption_spactrum.png' %fname)
        plt.close()
        #pl spectrum
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = 'intensity / a.u.')
        ax1.plot(PL_energy[0], PL_energy[1])
        plt.text(600, np.max(PL_energy[1]), 'absorbed_power : {:.3f} mW\nrelative_QE : {:.5f}'.format(uv_abs*1000, relative_QE))
        plt.xlim(300,800)
        plt.savefig('%s_PL_spactrum.png' %fname)
        plt.close()

    log('Wash flow cell')
    mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 1*ml, velocity=1000)
    for _ in range(4):
        mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', to_cell_over, velocity=500)

    #clean the flow cell
    log('Wash vial')
    for _ in range(3):
        mgr.pump.draw_and_dispense('valve', 'vial', 1.5*total_volume, velocity=1000)
        mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    return results


if __name__ == '__main__':


    # mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', 0.5*ml, velocity=500)
    # mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', 0.5*ml, velocity=500)
    
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', 1*ml, velocity=500)
    
    ##measurment#############
    import os 

    conc = 1
    sample = 'DSB'
    path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200805_test_abs_PL'
    if os.path.exists(path) == False:
        os.mkdir(path)
    fname = '%s/%suM_%s_in_ACN_6' %(path, conc, sample)

    measure_abs_PL2(fname, dilution = 2, dark_average = 50, abs_average=200, abs_exposure=0.015, led_power =5, PL_average=30, PL_exposure=1, filter_size=20, uv_average=10)
    #############


    ###washing###############

    #mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 1*ml, velocity=2000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_abs', 2*ml, velocity=1000)

    #mgr.pump.draw(volume = 0.5*ml, valve='ACN')
    #mgr.pump.set_velocity(100)
    #mgr.pump.dispense_all(valve='flow_cell_abs')
    ##############


    ###measure line#############

    # mgr.pump.draw_and_dispense('air', 'flow_cell_abs', 2*ml, velocity=1000)
    # mgr.pump.draw_and_dispense('flow_cell_abs', 'waste', 2*ml, velocity=1000)
  
    #measure_line('flow_cell_abs', fill='ACN')

    # mgr.pump.draw(volume = 0.5*ml, valve='ACN' )
    # mgr.pump.dispense(volume = 0.35*ml, valve='flow_cell_abs')
    # mgr.pump.dispense(volume = 0.15*ml, valve='flow_cell_abs')

    #################

    #measure_abs_PL_manual_injection(None, dark_average = 50, abs_average=200, abs_exposure=0.015, PL_average=30, PL_exposure=0.01, filter_size=20, uv_average=10)
    
    # measure_abs_PL_manual_injection(fname, abs_average=3000, abs_exposure=0.005, PL_average=30, PL_exposure=0.01, filter_size=20, uv_average=10)
    
    
    # measure_abs_PL_manual(None, 2, 'wash')
 #   vial_wash(1*ml)
#     mgr.pump.draw_and_dispense('vial', 'waste', 2*ml, velocity = 1000)
# #    mgr.pump.draw_and_dispense('valve', 'flow_cell', 2*ml, velocity = 500)

#     sample_info ={
#         'name' : 'di-nBu-PTPTP',
#         'concentration' : '2uM',
#         'solvent' : 'THF'
#     }

#     path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200714_test_autoPL/'
#     sample_name = '20uL-di-nBu-PTPTP_in_THF' 
#     measure_TE(path + sample_name, dilution = 1, **sample_info)

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






