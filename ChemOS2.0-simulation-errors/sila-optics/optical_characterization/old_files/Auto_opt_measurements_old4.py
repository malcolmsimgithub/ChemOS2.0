from pylab.manager import Logger
import numpy as np
import matplotlib.pyplot as plt
import os, glob, time, copy, math
from transient_emission import Transient_Emission
from absorption_and_PL import Abs_PL
import pickle
import json
import csv
from proc_image import ImageHandler
from pylab.instruments import PumpPSD8


filedir = os.path.dirname(os.path.realpath(__file__))
config_file = '%s/.config_Absorption_and_PL.dat' %filedir


## settings
ml = 1e-3

buffer_volume = 0.025*ml
loop_to_pump  = 0.04*ml # (under)

flow_cells = {
    'TE' : {
        'name' : 'flow_cell',
        'to_cell_under' : 0.46*ml,
        'to_cell_over' : 0.51*ml
    },
    'absorption' : {
        'name' : 'flow_cell_abs',
        'to_cell_under' : 0.245*ml,
        'to_cell_over' : 0.28*ml       
    },
    'PL' : {
        'name' : 'flow_cell_PL',
        'to_cell_under' : 0.34*ml,
        'to_cell_over' : 0.39*ml
    } 
}

for val in flow_cells.values():
    val['cell_volume'] = round(val['to_cell_over'] - val['to_cell_under'] + 2*buffer_volume, 6)


class Optical_measurements():

    def __init__(self, TE=True, AbsPL=True):

        self.log = Logger(stdout=True, logfile=None, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')

        pump_setting =   { 
                    "visa": "ASRL2::INSTR",
                    "addr": 0,
                    "syringe_volume": 1E-3,
                    "init_valve": 8,
                    "ports": {
                        "vial": 1,
                        "air": 2,
                        "flow_cell_PL": 3,
                        "flow_cell_abs" : 4,
                        "flow_cell": 5,
                        "water":6,
                        "waste": 7,
                        "valve": 8,
                    }
        }
        
        self.pump = PumpPSD8(**pump_setting)
        self.pump.set_velocity(1000)

        self.TE, self.AP = None, None
        if TE:           
            self.TE = Transient_Emission(device = True, DB = False)
            self.TE.detector_on()
        if AbsPL:
            self.abspl_config = self._load_abs_PL_config(config_file)
            self.AP = Abs_PL(led_power=self.abspl_config['PL']['led_power'], correction_datafile= self.abspl_config['data_processing']['correction_datafile'])
        self.IM = ImageHandler()


    def _load_abs_PL_config(self, config_file):
        
        with open(config_file) as content:
            config = json.loads(content.read())

        for val in config.values():
            for key in val.keys():
                if val[key] == 'True':
                    val[key] = True
                elif val[key] == 'False':
                    val[key] = False
                elif val[key] == 'None':
                    val[key] = None
        return config


    def _load_abs_PL_setting(self, **kwargs):
        import copy
        config = copy.deepcopy(self.abspl_config)
        setting = {}

        for val in config.values():
            for key in val.keys():
                if key in kwargs.keys():
                    val[key] = kwargs[key]
                setting[key] = val[key]

        return config, setting

    #overwrite pump methods
    def draw(self, valve, volume, velocity):
        self.pump.set_velocity(velocity)
        self.pump.draw(volume, valve = valve)

    def dispense(self, valve, volume, velocity):
        self.pump.set_velocity(velocity)
        self.pump.dispense(volume, valve = valve)

    def draw_full(self, valve, velocity):
        self.pump.set_velocity(velocity)
        self.pump.draw_full(valve = valve)

    def dispense_all(self, valve, velocity):
        self.pump.set_velocity(velocity)
        self.pump.dispense_all(valve=valve)

    def draw_and_dispense(self, from_valve, to_valve, volume, draw_velocity, dispense_velocity):
        self.draw(from_valve, volume, draw_velocity)
        self.dispense(to_valve, volume, dispense_velocity)


    def measure_line(self, line, fill=None):
        if fill:
            self.pump.draw_full(fill)
            action = self.pump.dispense
        else:
            action = self.pump.draw

        c = ''
        volume = 0
        while c != 'c':
            action(volume=0.01*ml, valve=line)
            volume += 0.01
            c = input(f'{volume:.3f}')

        self.pump.dispense_all('waste')


    def cell_wash(self, cells, volume, repeat):
        for cell in cells:
            self.log('Wash %s' %cell)
            for _ in range(repeat):
                self.draw_and_dispense('valve', cell, volume, draw_velocity=500, dispense_velocity=1000)


    def vial_wash(self, volume, repeat):
        self.log('Wash vial')
        for _ in range(repeat):
            self.draw_and_dispense('valve', 'vial', volume, draw_velocity = 1000, dispense_velocity = 1000)
            self.draw_and_dispense('vial', 'waste', 1*ml, draw_velocity = 2000, dispense_velocity = 2000)


    def fill_reference(self, water_ratio, measurements):

        num = len([measurement for measurement in measurements if measurement == 'absorption' or measurement == 'PL'])

        if water_ratio == 0 or num == 0:
            pass
        else:
            self.log('mixing reference solution')
            self.draw_and_dispense('valve', 'vial', 0.5 * num * (1-water_ratio)*ml, draw_velocity = 1000, dispense_velocity = 1000)
            self.draw_and_dispense('water', 'vial', 0.5 * num * water_ratio*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
            for _ in range(2):
                self.draw_and_dispense('vial', 'vial', 1*ml, draw_velocity = 2000, dispense_velocity = 2000)
            
            if 'absorption' in measurements:
                self.log('Sending reference solution to absorption cell')
                self.draw_and_dispense('vial', 'flow_cell_abs', 0.5*ml, draw_velocity = 1000, dispense_velocity = 1000)
            if 'PL' in measurements:
                self.log('Sending reference solution to PL cell')
                self.draw_and_dispense('vial', 'flow_cell_PL', 0.5*ml, draw_velocity = 1000, dispense_velocity = 1000)


    def sample_dilution(self, total_volume, dilution, setting):

        self.log('Discard valve to pump volumes')
        self.draw_and_dispense('valve', 'waste', loop_to_pump, setting['dilution_draw_velocity'], setting['dilution_dispense_velocity'])

        self.log(f'Draw {total_volume*1000:.1f} ml and mix in vial')
        self.draw_and_dispense('valve', 'vial', total_volume, setting['dilution_draw_velocity'], setting['dilution_dispense_velocity'])
        

        if dilution > 1:
            for _ in range(2):
                self.draw_and_dispense('vial', 'vial', total_volume *2, setting['dilution_draw_velocity'], setting['dilution_dispense_velocity']) # mix 3 times


    def send_sample_to_cell(self, cell, draw_velocity = 500, dispense_velocity = 1000):

        flow_cell = flow_cells[cell]
        volume = flow_cell['cell_volume']

        #make flow cell empty
        self.log('Emplying the flow cell {timestamp}')
        self.pump.draw_and_dispense('air', flow_cell['name'], 0.5*ml, velocity=1000) #better not to draw solvent to the pump
        self.pump.draw_and_dispense('air', flow_cell['name'], 1*ml, velocity=2000) #to prevent unintentioanl dilution
        #self.pump.draw_and_dispense(flow_cell['name'], 'waste', 1*ml, velocity=2000) #to prevent unintentioanl dilution
        
        #self.pump.draw_and_dispense(flow_cell['name'], 'waste', 2*ml, velocity=2000)

        self.log(f'Send {volume/ml:.3f} ml of solution to flow cell')
        #self.draw_full('vial', draw_velocity)
        draw_volume = min(1*ml, volume *3)
        self.draw('vial', draw_volume , draw_velocity)

        self.dispense('vial', draw_volume - volume, dispense_velocity)
        self.dispense_all(flow_cell['name'], dispense_velocity)

        #send sample to the flowcell
        self.draw_and_dispense('air', flow_cell['name'], flow_cell['to_cell_under']- buffer_volume, 500, dispense_velocity)


    def collect_sample_from_cell(self, cell, draw_velocity = 500, dispense_velocity = 1000):
        
        flow_cell = flow_cells[cell]
        volume = flow_cell['cell_volume']

        self.log(f'Collect {volume/ml:.3f} ml of solution to flow cell')
        self.draw(flow_cell['name'], flow_cell['to_cell_over'] + buffer_volume *2, draw_velocity)
        #self.draw_full(flow_cell['name'], draw_velocity)
        self.dispense_all('vial', dispense_velocity)



    def save_data(self, fname, data):
        np.savetxt(fname, data, delimiter=',', header='wl (nm), abs (%), pl (and abs/pl pair repeats)')


    def write_result_csv(self, fname, results):

        if results['absorption'] == None and results['PL'] == None:
            return

        header1 = []
        header2 = []
        stats = []
        d = []

        if results['absorption'] != None:
            header1.extend(['abs_max', 'abs_lambda_max'])
            stats.extend([results['absorption']['max'], results['absorption']['lambda_max'] ])
            if results['absorption']['end_wavelength']:
                header1.extend(['abs_end'])
                stats.extend([results['absorption']['end_wavelength']])
            header2.extend(['wavelength/nm', 'abs_ref', 'abs_sample', 'transmittance', 'absorbance'])
            d.extend([results['absorption']['reference'][0], 
                        results['absorption']['reference'][1], 
                        results['absorption']['sample'][1], 
                        results['absorption']['transmittance'][1], 
                        results['absorption']['absorbance'][1]])
            
        if results['PL'] != None:
            header1.extend(['uv_ref(W)', 'uv_absorption(W)', 'uv_absorption_maintenance','uv_absorption_maintenance(at_1min)', 'degradation_rate(s-1)', 'PL_max','PL_lambda_max', 'relative_QY'])
            stats.extend([results['uv']['reference'], 
                          results['uv']['absorption'], 
                          results['uv']['absorption_maintenance'],
                          results['uv']['absorption_maintenance(at_1min)'],
                          results['uv']['degradation_rate'],
                          results['PL']['max'], 
                          results['PL']['lambda_max'], 
                          results['PL']['relative_QY']])
            if 'max_gain_factor' in results['PL'].keys():
                header1.extend(['max_gain_factor(cm2 s)'])
                stats.extend([results['PL']['max_gain_factor']])
            header2.extend(['wavelength/nm', 'PL(energy/s/nm)', 'PL(photons/s/nm)', 'PL(photons/s/Hz)', 'gain_factor(cm2 s)'])
            d.extend([results['PL']['energy'][0], 
                        results['PL']['energy'][1], 
                        results['PL']['photons'][1],
                        results['PL']['freq_spectrum'][1],
                        results['PL']['gain_spectrum'][1]])

        with open(fname, 'w', newline = '') as f:

            writer = csv.writer(f)
            writer.writerow(header1)
            writer.writerow(stats)
            d = np.transpose(np.asarray(d))
            writer.writerow(header2)
            writer.writerows(d)
        
        print('%s was saved' %fname)



    def abspl_result_plots(self, fname, results):

        if results['absorption'] != None:
            text = 'abs_max : {:.4f} ({:.1f} nm)'.format(results['absorption']['max'], results['absorption']['lambda_max'])
            if results['absorption']['end_wavelength']:
                text = text + '\nabs_end : {:.1f} nm'.format(results['absorption']['end_wavelength'])
            self.AP.result_plot(results['absorption']['absorbance'], 'Absorbance', text=text, xrange = [300,800], show_plot=False, save_filename='%s_absorption_spactrum.png' %fname)
 
        if results['PL'] != None:
            text = 'PL_max : {:.4f} ({:.1f} nm)\nabsorbed_power : {:.3f} mW\nrelative_QY : {:.3f}'.format\
                           (results['PL']['max'],
                            results['PL']['lambda_max'],
                            results['uv']['absorption']*1000,
                            results['PL']['relative_QY'])
            if results['uv']['absorption_maintenance'] is not None:
                text = text + '\nuv_absorption_maintenance(%) : {:.2f}\nuv_absorption_maintenance_1min(%) : {:.2f}\ndegradation_rate(s-1) : {:.5f}'.format\
                            (results['uv']['absorption_maintenance']*100,
                            results['uv']['absorption_maintenance(at_1min)']*100,
                            results['uv']['degradation_rate'])
            self.AP.result_plot(results['PL']['energy'], 'Intensity/ (energy/s/nm)', text = text,\
                                                                xrange = [300,800], show_plot=False, save_filename='%s_PL_spactrum.png' %fname)
 
            self.AP.result_plot(results['PL']['gain_spectrum'], 'gain factor/ (cm2 s)', text = None,\
                                                                xrange = [300,800], show_plot=False, save_filename='%s_PL_gain_spactrum.png' %fname)

        if results['absorption'] != None and results['PL'] != None:
            self.AP.Abs_PL_plot(results['absorption']['absorbance'], results['PL']['energy'], text = None,\
                                                                xrange = [300,800], show_plot=False, save_filename='%s_Abs_PL_spactrum.png' %fname)

            if 'max_gain_factor' in results['PL'].keys():
                text = 'max_gain_factor(cm2 s) :\n {:.3e} ({:.1f} nm)'.format(results['PL']['max_gain_factor'], results['PL']['max_gain_wavelength'])
            else:
                text = None
            self.AP.Abs_PL_plot(results['absorption']['absorbance'], results['PL']['gain_spectrum'], text = text, xrange = [300,800], 
                                        ylabel = 'Nomalized Abs. and gain_factor.', show_plot=False, save_filename='%s_Abs_PL_gain_spectrum.png' %fname)


    def join_images(self, fname, filelist):
        if len(filelist) > 1:
            self.IM.tile_img(h_num = 3, save_filename = '%s.png' %fname, mergin = [0,0,0,0], del_files=True, filelist=filelist)


    def measure_PL(self, setting):

        uv, PL = {}, {}

        l_index, u_index = self.AP._to_index_range(self.AP.wl, setting['calc_range'])

        #measure reference
        self.log('Measure reference uv power {timestamp}')
        uv['reference'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

        self.send_sample_to_cell(cell = 'PL', draw_velocity=setting['PL_draw_velocity'], dispense_velocity=setting['PL_dispense_velocity'])

        # need sometime until the solution reaches equiliblium (bubble??)
        print('waiting %s s for equibliration' %setting['PL_equibliration_time'])
        time.sleep(setting['PL_equibliration_time']) 

        #measure PL spectrum
        self.log('Measure uv power {timestamp}')
        uv['before'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

        self.log('Adjust PL exposure {timestamp}')
        PL['exposure'] = self.AP.adjust_PL_exposure(setting['PL_initial_exposure'], setting['PL_max_exposure'], setting['PL_target_intensity'],\
                                                                        setting['led_power'], filter_size=setting['filter_size'], average = 5)
        self.log('Measure PL dark spectrum {timestamp}')
        self.AP.measure_dark_spectrum(20, PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
        self.log('Measure PL spectrum {timestamp}')
        start_time = time.time()
        PL_raw = self.AP.measure_PL_spectrum(setting['PL_average'],  PL['exposure'],  led_power = setting['led_power'],\
                                    filter_size = setting['filter_size'], dark_correction = True, spectral_correction = True, do_plot = False)
        end_time = time.time()
        self.log('Measure uv power {timestamp}')
        uv['after'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])


        #calculations 
        uv['absorption'] = uv['reference'] - (uv['before'] + uv['after'])/2
        PL['energy'] = np.asarray([PL_raw[0], PL_raw[1]/PL['exposure']])
        PL['photons'] = np.asarray([PL_raw[0], PL['energy'][1]*PL['energy'][0]*1E-9])
        PL['freq_spectrum'] = self.AP.calc_freq_spectrum(PL['photons'], calc_range= setting['calc_range'])
        PL['gain_spectrum'] = self.AP.calc_gain_spectrum(PL['photons'], calc_range= setting['calc_range'])
        PL['relative_QY'] = np.sum(PL['photons'][1][l_index:u_index+1])/uv['absorption']
        uv['absorption_maintenance'] = min(1, (uv['reference'] - uv['after']) /(uv['reference']-uv['before']))
        uv['degradation_rate'] = - math.log(uv['maintenance'])/(end_time-start_time)
        uv['absorption_maintenance(at_1min)'] = uv['absorption_maintenance'] * math.exp(-(60-(end_time-start_time))* uv['degradation_rate'])
        PL['max'], PL['lambda_max'] = self.AP.find_max(PL['energy'], analysis_range=setting['calc_range'])

        print('PL_max :  {:.4f} ({:.1f} nm)'.format(PL['max'], PL['lambda_max']))
        print('uv_ref : {:.3f} mW'.format(uv['reference'] * 1000))
        print('uv_sample_before : {:.3f} mW'.format(uv['before'] * 1000))
        print('uv_sample_after : {:.3f} mW'.format(uv['after'] * 1000))
        print('uv_absorption_maintenance : {:.3f}%'.format(uv['absorption_maintenance']*100))
        print('degradation_rate(s-1): {:.3f}'.format(uv['degradation_rate']))
        print('uv_absorption_maintenance(at_1min) : {:.3f}%'.format(uv['absorption_maintenance(at_1min)']*100))
        print('relative_QY : {:.3f}'.format(PL['relative_QY']))

        input('pause')
        self.collect_sample_from_cell(cell = 'PL', draw_velocity=setting['PL_draw_velocity'], dispense_velocity=setting['PL_dispense_velocity'])

        return uv, PL


    def measure_PL2(self, setting):

        uv, PL = {}, {}

        l_index, u_index = self.AP._to_index_range(self.AP.wl, setting['calc_range'])

        #measure reference
        self.log('Measure reference uv power {timestamp}')
        uv['reference'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

        self.send_sample_to_cell(cell = 'PL', draw_velocity=setting['PL_draw_velocity'], dispense_velocity=setting['PL_dispense_velocity'])

        # need sometime until the solution reaches equiliblium (bubble??)
        print('waiting %s s for equibliration' %setting['PL_equibliration_time'])
        time.sleep(setting['PL_equibliration_time']) 


        self.log('Adjust PL exposure {timestamp}')
        PL['exposure'] = self.AP.adjust_PL_exposure(setting['PL_initial_exposure'], setting['PL_max_exposure'], setting['PL_target_intensity'],\
                                                                        setting['led_power'], filter_size=setting['filter_size'], average = 5)
        self.log('Measure PL dark spectrum {timestamp}')
        self.AP.measure_dark_spectrum(20, PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
        self.log('Measure PL spectrum {timestamp}')
        results = self.AP.measure_PL_spectrum2(setting['PL_average'],  PL['exposure'],  led_power = setting['led_power'],\
                                    filter_size = setting['filter_size'], dark_correction = True, spectral_correction = True, do_plot = False)


        #calculations 
        uv['absorption'] = uv['reference'] - results['uv_average']
        PL['energy'] = np.asarray([results['PL_spectrum'][0], results['PL_spectrum'][1]/PL['exposure']])
        PL['photons'] = np.asarray([results['PL_spectrum'][0], PL['energy'][1]*PL['energy'][0]*1E-9])
        PL['freq_spectrum'] = self.AP.calc_freq_spectrum(PL['photons'], calc_range= setting['calc_range'])
        PL['gain_spectrum'] = self.AP.calc_gain_spectrum(PL['photons'], calc_range= setting['calc_range'])
        PL['relative_QY'] = np.sum(PL['photons'][1][l_index:u_index+1])/uv['absorption']
        PL['max'], PL['lambda_max'] = self.AP.find_max(PL['energy'], analysis_range=setting['calc_range'])

        
        print('PL_max :  {:.4f} ({:.1f} nm)'.format(PL['max'], PL['lambda_max']))
        print('uv_ref : {:.3f} mW'.format(uv['reference'] * 1000))
        print('uv_start : {:.3f} mW'.format(results['uv_start'] * 1000))
        print('uv_end : {:.3f} mW'.format(results['uv_end'] * 1000))
        print('relative_QY : {:.3f}'.format(PL['relative_QY']))

        if uv['absorption']/uv['reference'] > 0.003:
            uv['absorption_maintenance'] = min(1, (uv['reference']-results['uv_end']) /(uv['reference']-results['uv_start']))
            uv['degradation_rate'] = - math.log(uv['absorption_maintenance'])/results['duration']
            uv['absorption_maintenance(at_1min)'] = uv['absorption_maintenance'] * math.exp(-(60-results['duration'])* uv['degradation_rate'])

            print('uv_absorption_maintenance : {:.3f}%'.format(uv['absorption_maintenance']*100))
            print('degradation_rate: {:.3f}'.format(uv['degradation_rate']))
            print('uv_absorption_maintenance(at_1min) : {:.3f}%'.format(uv['absorption_maintenance(at_1min)']*100))
        else:
            uv['absorption_maintenance'] = None
            uv['degradation_rate'] = None
            uv['absorption_maintenance(at_1min)'] = None

        
        self.collect_sample_from_cell(cell = 'PL', draw_velocity=setting['PL_draw_velocity'], dispense_velocity=setting['PL_dispense_velocity'])

        return uv, PL


    def measure_absorption(self, setting):

        absorption = {}

        #measure reference
        self.log('Measure absorption reference {timestamp}')
        self.AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
        absorption['reference'] = self.AP.measure_transmission_spectrum(setting['abs_average'], setting['abs_exposure'], filter_size = setting['filter_size'], dark_correction=True, do_plot=False)
    

        self.send_sample_to_cell(cell = 'absorption', draw_velocity=setting['abs_draw_velocity'], dispense_velocity=setting['abs_dispense_velocity'])

        #need sometime until the solution reaches equiliblium (bubble??)
        print('waiting %s s for equibliration' %setting['abs_equibliration_time'])
        time.sleep(setting['abs_equibliration_time']) 

        #measure absorption spectrum
        self.log('Measure absorption dark {timestamp}')
        self.AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
        self.log('Measure absorption spectrum  {timestamp}')
        absorption['sample'], absorption['transmittance'], absorption['absorbance'] = self.AP.measure_absorption_spectrum(setting['abs_average'], setting['abs_exposure'], \
                                            absorption['reference'], filter_size=setting['filter_size'], dark_correction=True, do_plot=False)

        absorption['max'], absorption['lambda_max'] = self.AP.find_max(absorption['absorbance'], analysis_range=setting['calc_range']) 
        absorption['end_index'], absorption['end_wavelength'] = self.AP.find_abs_end(absorption['absorbance'], setting['absorption_threshold'], analysis_range=setting['calc_range'])     
        print('abs_max :  {:.4f} ({:.1f} nm)'.format(absorption['max'], absorption['lambda_max']))
        if absorption['end_wavelength']:
            print('abs_end :  {:.1f} nm'.format(absorption['end_wavelength']))

        input('pause')
        self.collect_sample_from_cell(cell = 'absorption', draw_velocity=setting['abs_draw_velocity'], dispense_velocity=setting['abs_dispense_velocity'])

        return absorption
    

        ## start pl and abs measurement with pump
    def measure_TE(self, fname, setting, **kwargs):

        self.send_sample_to_cell(cell = 'TE', draw_velocity=setting['TE_draw_velocity'], dispense_velocity=setting['TE_dispense_velocity'])

        self.log('Measure transient emission {timestamp}')
        #TE.measure_TE(**kwargs)

        result = self.TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)

        self.collect_sample_from_cell(cell = 'TE', draw_velocity=setting['TE_draw_velocity'], dispense_velocity=setting['TE_dispense_velocity'])
        #self.pump.draw_and_dispense('flow_cell', 'vial', 1*ml, velocity=1000)

        return result


    def _check_instruments(self, measurements):
        flg = True
        if 'absorption' in measurements and self.AP == None:
            print('absorption setup is not active.')
            flg = False    
        elif 'PL' in measurements and self.AP == None:
            print('PL setup is not active.')
            flg = False
        elif 'TE' in measurements and self.TE == None:
            print('TE setup is not active.')
            flg = False
        return flg


    def _find_maximum_gain(gain_spctrum , absorption_end_index, calc_range = [300,800]):

        l_index, u_index = self.AP._to_index_range(gain_spectrum[0], calc_range)
        spectrum = gain_spectrum[:, l_index:u_index+1]

        max_gain = np.max(spectrum[1][absorption_end_index : ])
        max_gain_wavelength = spectrum[0][np.argmax(spectrum[1][absorption_end_index : ])+ absorption_end_index]
        return max_gain, max_gain_wavelength


    def do_measurements(self, fname, measurements, dilution, water_ratio, **kwargs):

        if  self._check_instruments(measurements) == False:
            raise Exception()

        results, uv, PL, absorption, TE = {}, None, None, None, None
        
        config, setting = self._load_abs_PL_setting(**kwargs)

        # calculation to send to flow_cell
        total_volume = dilution*0.1*ml

        dilution_actual = total_volume/(0.1*ml)
        print('dilution : {:.2f}'.format(dilution_actual))

        self.fill_reference(water_ratio, measurements)

        self.sample_dilution(total_volume, dilution_actual, setting)

        time.sleep(1)

        #do measurements
        if 'absorption' in measurements:
            absorption = self.measure_absorption(setting)

        #TODO dilution for PL

        if 'PL' in measurements:
            uv, PL = self.measure_PL2(setting)
            
        if 'TE' in measurements:
            try:
                TE = self.measure_TE(fname, setting, **kwargs)
                os.remove('%s_TE.pkl' %fname)
                del(TE['metadata']['sample'])
                config['TE'].update(copy.deepcopy(TE['metadata']))
                del(TE['metadata'])
            except:
                print('TE measurement failed')
            
        if absorption is not None and 'gain_spectrum' in PL.keys():
            if absorption['end_index']:
                PL['max_gain_factor'], PL['max_gain_wavelength'] = self._find_maximum_gain(PL['gain_spectrum'], absorption['end_index'], calc_range = setting['calc_range'])

        results = {'uv' : uv,  'absorption' : absorption, 'PL' :  PL,  'TE' : TE, 'metadata' : config}

        if fname:
            self.AP.save_pkl(fname + '.pkl', results)
            self.write_result_csv(fname + '_AbsPL.csv', results)
            self.abspl_result_plots(fname, results)

        self.join_images(fname, glob.glob('%s*.png' %fname))

        self.log('discard sample {timestamp}')
        self.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)

        #wash flow cells and vial
        for i, cell in enumerate(measurements):
            if i == 0:
                n = 3
            else: n = 2
            self.cell_wash([flow_cells[cell]['name']], 0.5*ml, n)

        self.vial_wash(1.5*total_volume, 4)

        return results


if __name__ == '__main__':


    #initialize the device############
    Opt = Optical_measurements(TE = True, AbsPL = True)
    ##############

    #measurment_abs_PL_TE#############
    import os 

    sample_info ={
        'name' : 'DSB',
        'concentration(uM)' : 1,
        'solvent' : 'ACN'
    }

    path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200817_test_abs_PL_TE_2'
    if os.path.exists(path) == False:
       os.mkdir(path)
    fname = '%s/%s_%s_in_%s_3' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])

    Opt.do_measurements(fname, measurements = ['absorption', 'PL','TE'], dilution = 1, water_ratio = 0, **sample_info)
    #Opt.do_measurements(fname, measurements = ['TE'], dilution = 1, water_ratio = 0, **sample_info)
    ############

    #fill cells########
    #Opt.draw_and_dispense('valve', 'flow_cell_abs', 1*ml, draw_velocity=500, dispense_velocity=500)
    #Opt.draw_and_dispense('valve', 'flow_cell_PL', 1*ml, draw_velocity=500, dispense_velocity=500)
    #Opt.draw_and_dispense('valve', 'flow_cell', 1*ml, draw_velocity=500, dispense_velocity=500)
    #Opt.vial_wash(0.15*ml, 2)
    ############

     

    ###measure line#############
    #Opt = Optical_measurements(TE=False, AbsPL=True)

    #Opt.mgr.pump.draw_and_dispense('air', 'flow_cell_PL', 1*ml, velocity=1000)
    #Opt.mgr.pump.draw_and_dispense('flow_cell_PL', 'waste', 1*ml, velocity=1000)
    #Opt.measure_line('flow_cell_PL', fill='valve')
    #Opt.mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.34*ml, velocity=1000)
    #input()
    #Opt.mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.39*ml, velocity=1000)
    
    ###########################

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


    ###test new PL#############
    # results = Opt.AP.measure_PL_spectrum2(30,  0.3,  led_power = 10,\
    #                                 filter_size = 20, dark_correction = False, spectral_correction = True, do_plot = False)
    # print(results['uv_average'], results['uv_start'], results['uv_end'], results['duration'], results['average'] )
    # plt.plot(results['PL_spectrum'][0], results['PL_spectrum'][1])
    # plt.show()
    ##############





