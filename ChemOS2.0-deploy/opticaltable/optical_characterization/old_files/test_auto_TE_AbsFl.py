from pylab.manager import Manager, Logger
import pylab.instruments as ins
import numpy as np
import matplotlib.pyplot as plt
import os, time
from transient_emission import Transient_Emission
from absorption_and_PL import Abs_PL
import pickle
import json
import csv


filedir = os.path.dirname(os.path.realpath(__file__))
config_file = '%s/.config_Absorption_and_PL.dat' %filedir

## settings
ml = 1e-3
loop_to_pump  = 0.04*ml # (under)
to_cell_under = 0.46*ml #0.43
to_cell_over  = 0.51*ml #0.57

to_cell_under_PL = 0.35*ml 
to_cell_over_PL  = 0.50*ml
# to_cell_under_PL = 0.46*ml 
# to_cell_over_PL  = 0.63*ml


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
                "flow_cell_abs" : 2,
                "air": 3,
                "flow_cell_PL": 4,
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


def vial_wash(volume, repeat):
    for _ in range(repeat):
        mgr.pump.draw_and_dispense('valve', 'vial', volume, velocity=1000)
        mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)

def cell_wash(cells, volume, repeat):
    for cell in cells:
        for _ in range(repeat):
            mgr.pump.draw_and_dispense('valve', cell, volume, velocity=500)

## start pl and abs measurement with pump
def measure_TE(fname, dilution, **kwargs):

    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    # mgr.valve.move_to('A')
    # _ = input('Please load your material')

    # # switch position
    # log('Switch valve to inject position')
    # mgr.valve.move_to('B')

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



def _load_abs_PL_setting(config_file, **kwargs):
    
    with open(config_file) as content:
        config = json.loads(content.read())

    setting = {}

    for val in config.values():
        for key in val.keys():
            if key in kwargs.keys():
                val[key] = kwargs[key]
            if val[key] == 'True':
                val[key] = True
            elif val[key] == 'False':
                val[key] = False
            elif val[key] == 'None':
                val[key] = None
            setting[key] = val[key]
    return config, setting


def write_result_csv(fname, results):
    with open(fname, 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(['uv_ref(W)', 'uv_sample_before(W)', 'uv_sample_after(W)', 'uv_change(%)', 'uv_absorption(W)',\
                                             'abs_max', 'abs_lambda_max', 'PL_max','PL_lambda_max', 'relative_QY'])
        writer.writerow([results['uv']['reference'], 
                         results['uv']['before'], 
                         results['uv']['after'], 
                         results['uv']['change'], 
                         results['uv']['absorption'], 
                         results['absoption']['max'],
                         results['absoption']['lambda_max'],
                         results['PL']['max'],
                         results['PL']['lambda_max'],
                         results['PL']['relative_QY']])
        d = np.transpose(np.asarray([results['absorption']['reference'][0], 
                                     results['absorption']['reference'][1], 
                                     results['absorption']['sample'][1], 
                                     results['absorption']['transmittance'][1], 
                                     results['absorption']['absorbance'][1], 
                                     results['PL']['energy'][1], 
                                     results['PL']['photons'][1]]))
        writer.writerow(['wavelength/nm', 'abs_ref', 'abs_sample', 'transmittance', 'absorbance', 'PL(energy/s/nm)', 'PL(photons/s/nm)'])
        writer.writerows(d)
    print('%s was saved' %fname)


def result_plots(fname, results):

    #absoption spectrum
    AP.result_plot(results['absorption']['absorbance'], 'Absorbance', xrange = [300,800], show_plot=False, save_filename='%s_absorption_spactrum.png' %fname)
    #pl spectrum
    text = 'absorbed_power : {:.3f} mW\nrelative_QY : {:.3f}\nuv_change(%) : {:.2f}'.format(results['uv']['absorption']*1000, \
                                                                        results['PL']['relative_QY']*1000, results['uv']['change']*100)
    AP.result_plot(results['PL']['energy'], 'Intensity/ (energy/s/nm)', text = text,\
                                                        xrange = [300,800], show_plot=False, save_filename='%s_PL_spactrum.png' %fname)

    AP.Abs_PL_plot(results['absorption']['absorbance'], results['PL']['energy'], text = None,\
                                                        xrange = [300,800], show_plot=False, save_filename='%s_Abs_PL_spactrum.png' %fname)



def show_results(results):

    print('abs_max :  {:.4f} ({:.1f} nm)'.format(results['absorption']['max'], results['absorption']['lambda_max']))
    print('PL_max :  {:.4f} ({:.1f} nm)'.format(results['PL']['max'], results['PL']['lambda_max']))
    print('uv_ref : {:.3f} mW'.format(results['uv']['reference'] * 1000))
    print('uv_sample : {:.3f} mW'.format(results['uv']['sample'] * 1000))
    print('uv_change : {:.3f}%'.format(results['uv']['change']))
    print('relative_QY : {:.3f}'.format(results['PL']['relative_QY']))



def measure_abs_PL(fname, dilution, **kwargs):
    #currently cell volume is 140uL so dilution should be more than 2

    results, uv, PL, absorption = {}, {}, {}, {}

    # calculation to send to flow_cell
    buffer_volume = 0.025*ml
    flow_cell_volume = to_cell_over_PL - to_cell_under_PL+ 2*buffer_volume
    total_volume = max(dilution*0.1*ml, flow_cell_volume)

    dilution_actual = total_volume/(0.1*ml)

    print('dilution : {:.2f}'.format(dilution_actual))

    #initialize the measurement device
    config, setting = _load_abs_PL_setting(config_file, **kwargs)
    AP = Abs_PL(led_power=setting['led_power'], correction_datafile= setting['correction_datafile'])

    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=100)

    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)

    if dilution > 1:
        for i in range(2):
            mgr.pump.draw_and_dispense('vial', 'vial', total_volume *2, velocity =1000) # mix 3 times

    #measure backgrounds
    time.sleep(5)  # need sometime until the solution reaches equiliblium

    log('Measure uv power {timestamp}')
    uv['reference'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

    log('Measure background spectrum {timestamp}')
    AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
    absorption['reference'] = AP.measure_transmission_spectrum(setting['abs_average'], setting['abs_exposure'], filter_size = setting['filter_size'], dark_correction=True, do_plot=False)
    #input('pause')


    #draw solvent from the flow cell
    log('Discard solvent {timestamp}')
    mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 2*ml, velocity=2000)


    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')

    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='waste')
    mgr.pump.set_velocity(1000)
    mgr.pump.dispense_all(valve='flow_cell_PL')

    #send sample to the flowcell
    mgr.pump.draw(volume = to_cell_under_PL- buffer_volume, valve = 'air')
    mgr.pump.dispense_all(valve = 'flow_cell_PL')

    # need sometime until the solution reaches equiliblium
    print('waiting %s s for equibliration' %setting['equibliration_time'])
    time.sleep(setting['equibliration_time']) 

    #AP.equilibration(abs_average, abs_exposure,  equilibration_time = 10, threshold = 0.002, wavelength = 330, time_step = 1, time_out= 60, filter_size = filter_size)

    #measurement
    log('Measure absorbance {timestamp}')
    AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
    absorption['sample'], absorption['transimittance'], absorption['absorbance'] = AP.measure_absorption_spectrum(setting['abs_average'], setting['abs_exposure'], \
                                        absorption['reference'], filter_size=setting['filter_size'], dark_correction=True, do_plot=False)


    log('Measure uv power {timestamp}')
    uv['before'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

    log('Adjust PL exposure {timestamp}')
    PL['exposure'] = AP.adjust_PL_exposure(setting['PL_initial_exposure'], setting['PL_max_exposure'], setting['PL_target_intensity'],\
                                                                    setting['led_power'], filter_size=setting['filter_size'], average = 5)
    log('Measure PL {timestamp}')
    AP.measure_dark_spectrum(5, PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
    #AP.measure_dark_spectrum(setting['PL_average'], PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
    PL_raw = AP.measure_PL_spectrum(setting['PL_average'],  PL['exposure'],  led_power = setting['led_power'],\
                                 filter_size = setting['filter_size'], dark_correction = True, spectral_correction = True, do_plot = False)

    log('Measure uv power {timestamp}')
    uv['after'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])


    #calculations
    uv['absorption'] = uv['reference'] - (uv['before'] + uv['after'])/2
    PL['energy'] = np.asarray([PL_raw[0], PL_raw[1]/PL['exposure']])
    PL['photons'] = np.asarray([PL_raw[0], PL['energy'][1]*PL['energy'][0]*1E-9])
    PL['relative_QY'] = np.sum(PL['photons'][1])/uv['absorption']
    uv['change'] = 100 * (uv['after']-uv['before'])/uv['before']

    absorption['max'], absorption['lambda_max'] = AP.find_max(absorption, analysis_range=[300,800])
    PL['max'], PL['lambda_max'] = AP.find_max(PL['energy'], analysis_range=[300,800])


    results = {'uv' : uv,  'absorption' : absorption, 'PL' :  PL,  'metadata' : config}
    show_results(results)
    
    if fname:
        AP.save_pkl(fname + '.pkl', results)
        write_result_csv(fname + '.csv', results)
        result_plots(fname, results)

    log('Wash flow cell')
    mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 1*ml, velocity=1000)
    cell_wash(['flow_cell_PL'], to_cell_over_PL, 4)
    log('Wash vial')
    vial_wash(1.5*total_volume, 3)

    mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    return results



def measure_abs_PL_TE(fname, dilution, **kwargs):
    #currently cell volume is 140uL so dilution should be more than 2

    results, uv, PL, absorption = {}, {}, {}, {}
    
    # calculation to send to flow_cell
    buffer_volume = 0.025*ml
    flow_cell_volume = to_cell_over_PL - to_cell_under_PL+ 2*buffer_volume
    total_volume = max(dilution*0.1*ml, flow_cell_volume)
    dilution_actual = total_volume/(0.1*ml)

    print('dilution : {:.2f}'.format(dilution_actual))

    #initialize the measurement device
    TE = Transient_Emission(device = True, DB = False)
    TE.detector_on()

    config, setting = _load_abs_PL_setting(config_file, **kwargs)
    AP = Abs_PL(led_power=setting['led_power'], correction_datafile= setting['correction_datafile'])


    log('Discard valve to pump volumes')
    mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=100)

    log(f'Draw {dilution*0.1:.1f} ml and mix in vial')
    mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)

    if dilution > 1:
        for i in range(2):
            mgr.pump.draw_and_dispense('vial', 'vial', total_volume *2, velocity =1000) # mix 3 times

    #measure backgrounds
    time.sleep(5)  # need sometime until the solution reaches equiliblium

    log('Measure uv power {timestamp}')
    uv['reference'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

    log('Measure background spectrum {timestamp}')
    AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
    absorption['reference'] = AP.measure_transmission_spectrum(setting['abs_average'], setting['abs_exposure'], filter_size = setting['filter_size'], dark_correction=True, do_plot=False)
    #input('pause')


    #draw solvent from the flow cell
    log('Discard solvent {timestamp}')
    mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 2*ml, velocity=2000)


    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    mgr.pump.draw_full(valve='vial')
    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='vial')
    mgr.pump.set_velocity(1000)
    mgr.pump.dispense_all(valve='flow_cell_PL')

    #send sample to the flowcell
    mgr.pump.draw(volume = to_cell_under_PL- buffer_volume, valve = 'air')
    mgr.pump.dispense_all(valve = 'flow_cell_PL')

    # need sometime until the solution reaches equiliblium
    print('waiting %s s for equibliration' %setting['equibliration_time'])
    time.sleep(setting['equibliration_time']) 


    #measurement
    log('Measure absorbance {timestamp}')
    AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
    absorption['sample'], absorption['transimittance'], absorption['absorbance'] = AP.measure_absorption_spectrum(setting['abs_average'], setting['abs_exposure'], \
                                        absorption['reference'], filter_size=setting['filter_size'], dark_correction=True, do_plot=False)


    log('Measure uv power {timestamp}')
    uv['before'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

    log('Adjust PL exposure {timestamp}')
    PL['exposure'] = AP.adjust_PL_exposure(setting['PL_initial_exposure'], setting['PL_max_exposure'], setting['PL_target_intensity'],\
                                                                    setting['led_power'], filter_size=setting['filter_size'], average = 5)
    log('Measure PL {timestamp}')
    AP.measure_dark_spectrum(5, PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
    #AP.measure_dark_spectrum(setting['PL_average'], PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
    PL_raw = AP.measure_PL_spectrum(setting['PL_average'],  PL['exposure'],  led_power = setting['led_power'],\
                                 filter_size = setting['filter_size'], dark_correction = True, spectral_correction = True, do_plot = False)

    log('Measure uv power {timestamp}')
    uv['after'] = AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])


    #calculations
    uv['absorption'] = uv['reference'] - (uv['before'] + uv['after'])/2
    PL['energy'] = np.asarray([PL_raw[0], PL_raw[1]/PL['exposure']])
    PL['photons'] = np.asarray([PL_raw[0], PL['energy'][1]*PL['energy'][0]*1E-9])
    PL['relative_QY'] = np.sum(PL['photons'][1])/uv['absorption']
    uv['change'] = 100 * (uv['after']-uv['before'])/uv['before']

    absorption['max'], absorption['lambda_max'] = AP.find_max(absorption, analysis_range=[300,800])
    PL['max'], PL['lambda_max'] = AP.find_max(PL['energy'], analysis_range=[300,800])


    results = {'uv' : uv,  'absorption' : absorption, 'PL' :  PL,  'metadata' : config}
    show_results(results)
    
    if fname:
        AP.save_pkl(fname + '.pkl', results)
        write_result_csv(fname + '.csv', results)
        result_plots(fname, results)
        

    ###TE measurement#############
    mgr.pump.draw_full(valve='flow_cell_PL')

    flow_cell_volume = to_cell_over - to_cell_under + 2*buffer_volume
    log(f'Send {flow_cell_volume/ml:.2f} ml of solution to flow cell')
    #mgr.pump.draw_full(valve='vial')
    mgr.pump.dispense(volume=1*ml - flow_cell_volume, valve='vial')
    mgr.pump.dispense_all(valve='flow_cell')
    #send sample to the flowcell
    mgr.pump.draw(volume = to_cell_under - buffer_volume, valve = 'air')
    mgr.pump.dispense_all(valve = 'flow_cell')

    #input('pause')

    log('Measure transient emission {timestamp}')
    #TE.measure_TE(**kwargs)

    result = TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)


    log('Wash flow cells')
    mgr.pump.draw_and_dispense('flow_cell', 'waste', 1*ml, velocity=2000)
    mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)

    cell_wash(['flow_cell_PL', 'flow_cell'], 0.5*ml, 4)

    #clean the flow cell
    log('Wash vial')
    vial_wash(1.5*total_volume, 3)

    mgr.pump.draw_and_dispense('valve', 'waste', 0.5*ml, velocity=500)

    return results


if __name__ == '__main__':


    #mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.5*ml, velocity=500)
    #mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.5*ml, velocity=500)
    
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_PL', 1*ml, velocity=500)
    
    ##measurment_abs_PL#############
    # import os 

    # conc = 50
    # sample = 'DSB'
    # path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200807_test_abs_PL'
    # if os.path.exists(path) == False:
    #    os.mkdir(path)
    # fname = '%s/%suM_%s_in_ACN_3' %(path, conc, sample)

    # measure_abs_PL(fname, dilution = 2, equibliration_time = 60)
    #############

    #measurment_abs_PL_TE#############
    import os 

    sample_info ={
        'name' : 'diphenylantracene',
        'concentration' : '20uM',
        'solvent' : 'ACN'
    }

    path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200807_test_abs_PL_TE'
    if os.path.exists(path) == False:
       os.mkdir(path)
    fname = '%s/%s_%s_in_%s_4' %(path, sample_info['concentration'], sample_info['name'], sample_info['solvent'])

    measure_abs_PL_TE(fname, dilution = 2, **sample_info)
    ############

    # ###measurement TE#############
    # sample_info ={
    #     'name' : 'DSB',
    #     'concentration' : '10uM',
    #     'solvent' : 'ACN'
    # }
    # path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200807_test_abs_PL_TE'
    # if os.path.exists(path) == False:
    #    os.mkdir(path) 
    # fname = '%s/%s_%s' %(path, sample_info['name'], sample_info['concentration'])
    # measure_TE(fname, dilution = 2, **sample_info)
    ############

    ###washing##############

    # for _ in range(4):
    #     mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', to_cell_over_PL, velocity=500)
    # log('Wash Absorption flow cell')
    # for _ in range(4):
    #     mgr.pump.draw_and_dispense('valve', 'flow_cell', to_cell_over, velocity=500)


    ###washing###############

    #mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 1*ml, velocity=2000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_PL', 2*ml, velocity=1000)

    #mgr.pump.draw(volume = 0.5*ml, valve='ACN')
    #mgr.pump.set_velocity(100)
    #mgr.pump.dispense_all(valve='flow_cell_PL')
    ##############


    ###measure line#############

    # mgr.pump.draw_and_dispense('air', 'flow_cell_PL', 2*ml, velocity=1000)
    # mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 2*ml, velocity=1000)
  
    #measure_line('flow_cell_PL', fill='ACN')

    # mgr.pump.draw(volume = 0.5*ml, valve='ACN' )
    # mgr.pump.dispense(volume = 0.35*ml, valve='flow_cell_PL')
    # mgr.pump.dispense(volume = 0.15*ml, valve='flow_cell_PL')

    #################

    #measure_abs_PL_manual_injection(None, dark_average = 50, abs_average=200, abs_exposure=0.015, PL_average=30, PL_exposure=0.01, filter_size=20, uv_average=10)
    
    # measure_abs_PL_manual_injection(fname, abs_average=3000, abs_exposure=0.005, PL_average=30, PL_exposure=0.01, filter_size=20, uv_average=10)
    
    
    # measure_abs_PL_manual(None, 2, 'wash')
 #   vial_wash(1*ml)
#     mgr.pump.draw_and_dispense('vial', 'waste', 2*ml, velocity = 1000)
# #    mgr.pump.draw_and_dispense('valve', 'flow_cell', 2*ml, velocity = 500)



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

    #measure_line('flow_cell_PL', fill='ACN')
    #mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 2*ml, velocity = 1000)
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_PL', 0.63*ml, velocity = 1000)

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






