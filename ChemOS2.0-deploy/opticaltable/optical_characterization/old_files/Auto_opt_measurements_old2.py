from pylab.manager import Manager, Logger
import numpy as np
import matplotlib.pyplot as plt
import os, glob, time, copy
from transient_emission import Transient_Emission
from absorption_and_PL import Abs_PL
import pickle
import json
import csv
from proc_image import ImageHandler
#from pylab.instruments import PumpPSD8


filedir = os.path.dirname(os.path.realpath(__file__))
config_file = '%s/.config_Absorption_and_PL.dat' %filedir


## settings
ml = 1e-3

buffer_volume = 0.025*ml
loop_to_pump  = 0.04*ml # (under)

to_cell_under = 0.46*ml #0.43
to_cell_over  = 0.51*ml #0.57
flow_cell_volume = to_cell_over - to_cell_under+ 2*buffer_volume

to_cell_under_PL = 0.35*ml 
to_cell_over_PL  = 0.50*ml
flow_cell_volume_PL = to_cell_over_PL - to_cell_under_PL+ 2*buffer_volume


to_cell_under_abs = 0.235*ml 
to_cell_over_abs  = 0.270*ml
flow_cell_volume_abs = to_cell_over_abs - to_cell_under_abs+ 2*buffer_volume

flow_cells = {
    'TE' : {
        'name' : 'flow_cell',
        'to_cell_under' : 0.46*ml,
        'to_cell_over' : 0.51*ml
    },
    'absorption' : {
        'name' : 'flow_cell_abs',
        'to_cell_under' : 0.235*ml,
        'to_cell_over' : 0.27*ml       
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
        self.mgr = Manager([

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
            }
        ])

        self.mgr.pump.set_velocity(1000)
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

    def measure_line(self, line, fill=None):
        if fill:
            self.mgr.pump.draw_full(fill)
            action = self.mgr.pump.dispense
        else:
            action = self.mgr.pump.draw

        c = ''
        volume = 0
        while c != 'c':
            action(volume=0.01*ml, valve=line)
            volume += 0.01
            c = input(f'{volume:.3f}')

        self.mgr.pump.dispense_all('waste')


    def save_data(self, fname, data):
        np.savetxt(fname, data, delimiter=',', header='wl (nm), abs (%), pl (and abs/pl pair repeats)')


    def vial_wash(self, volume, repeat):
        for _ in range(repeat):
            self.mgr.pump.draw_and_dispense('valve', 'vial', volume, velocity=1000)
            self.mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)


    def cell_wash(self, cells, volume, repeat):
        for cell in cells:
            for _ in range(repeat):
                self.mgr.pump.draw_and_dispense('valve', cell, volume, velocity=500)


    def sample_dilution(self, total_volume, dilution):

        self.log('Discard valve to pump volumes')
        self.mgr.pump.draw_and_dispense('valve', 'waste', loop_to_pump, velocity=100)

        self.log(f'Draw {total_volume*1000:.1f} ml and mix in vial')
        self.mgr.pump.draw_and_dispense('valve', 'vial', total_volume, velocity=500)
        
        #self.mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', total_volume, velocity=500)
        #print("waiting 30s")
        #time.sleep(30)

        if dilution > 1:
            for _ in range(2):
                self.mgr.pump.draw_and_dispense('vial', 'vial', total_volume *2, velocity =1000) # mix 3 times


    def send_sample_to_cell(self, cell):

        flow_cell = flow_cells[cell]
        volume = flow_cell['cell_volume']

        #make flow cell empty
        self.log('Discard solvent {timestamp}')
        self.mgr.pump.draw_and_dispense('air', flow_cell['name'], 1*ml, velocity=1000) #better not to draw solvent to the pump
        self.mgr.pump.draw_and_dispense('air', flow_cell['name'], 1*ml, velocity=2000) #to prevent unintentioanl dilution
        #self.mgr.pump.draw_and_dispense(flow_cell['name'], 'waste', 1*ml, velocity=2000) #to prevent unintentioanl dilution
        
        #self.mgr.pump.draw_and_dispense(flow_cell['name'], 'waste', 2*ml, velocity=2000)


        self.log(f'Send {volume/ml:.3f} ml of solution to flow cell')
        self.mgr.pump.draw_full(valve='vial')
        self.mgr.pump.dispense(volume=1*ml - volume, valve='vial')
        self.mgr.pump.set_velocity(1000)
        self.mgr.pump.dispense_all(valve=flow_cell['name'])

        #send sample to the flowcell
        self.mgr.pump.draw(volume = flow_cell['to_cell_under']- buffer_volume, valve = 'air')
        #self.mgr.pump.draw(volume = flow_cell['to_cell_under']- buffer_volume, valve = 'valve')
        self.mgr.pump.dispense_all(valve = flow_cell['name'])



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
            header2.extend(['wavelength/nm', 'abs_ref', 'abs_sample', 'transmittance', 'absorbance'])
            d.extend([results['absorption']['reference'][0], 
                        results['absorption']['reference'][1], 
                        results['absorption']['sample'][1], 
                        results['absorption']['transmittance'][1], 
                        results['absorption']['absorbance'][1]])
            
        if results['PL'] != None:
            header1.extend(['uv_ref(W)', 'uv_sample_before(W)', 'uv_sample_after(W)', 'uv_change(%)', 'uv_absorption(W)','PL_max','PL_lambda_max', 'relative_QY'])
            stats.extend([results['uv']['reference'], results['uv']['before'], results['uv']['after'], results['uv']['change'], results['uv']['absorption'],\
                                                                                    results['PL']['max'], results['PL']['lambda_max'], results['PL']['relative_QY']])
            header2.extend(['wavelength/nm', 'PL(energy/s/nm)', 'PL(photons/s/nm)'])
            d.extend([results['PL']['energy'][0], 
                        results['PL']['energy'][1], 
                        results['PL']['photons'][1]])

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
            self.AP.result_plot(results['absorption']['absorbance'], 'Absorbance', text=text, xrange = [300,800], show_plot=False, save_filename='%s_absorption_spactrum.png' %fname)
        if results['PL'] != None:
            text = 'PL_max : {:.4f} ({:.1f} nm)\nabsorbed_power : {:.3f} mW\nrelative_QY : {:.3f}\nuv_change(%) : {:.2f}'.format\
                           (results['PL']['max'], results['PL']['lambda_max'], results['uv']['absorption']*1000, results['PL']['relative_QY'], results['uv']['change'])
            self.AP.result_plot(results['PL']['energy'], 'Intensity/ (energy/s/nm)', text = text,\
                                                                xrange = [300,800], show_plot=False, save_filename='%s_PL_spactrum.png' %fname)
        if results['absorption'] != None and results['PL'] != None:
            self.AP.Abs_PL_plot(results['absorption']['absorbance'], results['PL']['energy'], text = None,\
                                                                xrange = [300,800], show_plot=False, save_filename='%s_Abs_PL_spactrum.png' %fname)



    def show_abspl_results(self, results):

        if results['absorption'] != None:
            print('abs_max :  {:.4f} ({:.1f} nm)'.format(results['absorption']['max'], results['absorption']['lambda_max']))
        if results['uv'] != None:
            print('PL_max :  {:.4f} ({:.1f} nm)'.format(results['PL']['max'], results['PL']['lambda_max']))
            print('uv_ref : {:.3f} mW'.format(results['uv']['reference'] * 1000))
            print('uv_sample_before : {:.3f} mW'.format(results['uv']['before'] * 1000))
            print('uv_sample_after : {:.3f} mW'.format(results['uv']['after'] * 1000))
            print('uv_change : {:.3f}%'.format(results['uv']['change']))
            print('relative_QY : {:.3f}'.format(results['PL']['relative_QY']))


    def join_images(self, fname, filelist):
        if len(filelist) > 1:
            self.IM.tile_img(h_num = 2, save_filename = '%s.png' %fname, mergin = [0,0,0,0], del_files=True, filelist=filelist)


    def measure_PL(self, setting):

        uv, PL = {}, {}

        l_index, u_index = self.AP._to_index_range(self.AP.wl, setting['calc_range'])

        #measure reference
        self.log('Measure uv power {timestamp}')
        uv['reference'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

        self.send_sample_to_cell(cell = 'PL')

        # need sometime until the solution reaches equiliblium
        print('waiting %s s for equibliration' %setting['PL_equibliration_time'])
        time.sleep(setting['PL_equibliration_time']) 

        #measure PL spectrum
        self.log('Measure uv power {timestamp}')
        uv['before'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])

        self.log('Adjust PL exposure {timestamp}')
        PL['exposure'] = self.AP.adjust_PL_exposure(setting['PL_initial_exposure'], setting['PL_max_exposure'], setting['PL_target_intensity'],\
                                                                        setting['led_power'], filter_size=setting['filter_size'], average = 5)
        self.log('Measure PL {timestamp}')
        self.AP.measure_dark_spectrum(20, PL['exposure'], filter_size=setting['filter_size'],  do_plot=False) #dark for PL
        PL_raw = self.AP.measure_PL_spectrum(setting['PL_average'],  PL['exposure'],  led_power = setting['led_power'],\
                                    filter_size = setting['filter_size'], dark_correction = True, spectral_correction = True, do_plot = False)

        self.log('Measure uv power {timestamp}')
        uv['after'] = self.AP.measure_uv_power(counts = setting['uv_average'], led_power = setting['led_power'])


        #calculations
        uv['absorption'] = uv['reference'] - (uv['before'] + uv['after'])/2
        PL['energy'] = np.asarray([PL_raw[0], PL_raw[1]/PL['exposure']])
        PL['photons'] = np.asarray([PL_raw[0], PL['energy'][1]*PL['energy'][0]*1E-9])
        PL['relative_QY'] = np.sum(PL['photons'][1][l_index:u_index+1])/uv['absorption']
        uv['change'] = 100 * (uv['after'] -uv['before']) /(uv['reference']-uv['before'])
        PL['max'], PL['lambda_max'] = self.AP.find_max(PL['energy'], analysis_range=[300,800])

        #input('pause')

        self.mgr.pump.draw_and_dispense('flow_cell_PL', 'vial', 1*ml, velocity=1000)

        return uv, PL


    def measure_absorption(self, setting):

        absorption = {}

        #measure reference
        self.log('Measure absorption reference {timestamp}')
        self.AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
        absorption['reference'] = self.AP.measure_transmission_spectrum(setting['abs_average'], setting['abs_exposure'], filter_size = setting['filter_size'], dark_correction=True, do_plot=False)
        
#        input('pause')  #for manual measurements

        self.send_sample_to_cell(cell = 'absorption')

        #need sometime until the solution reaches equiliblium
        print('waiting %s s for equibliration' %setting['abs_equibliration_time'])
        time.sleep(setting['abs_equibliration_time']) 

        #measure absorption spectrum
        self.log('Measure absorbance {timestamp}')
        self.AP.measure_dark_spectrum(setting['dark_average'], setting['abs_exposure'], filter_size=setting['filter_size'],  do_plot=False)
        absorption['sample'], absorption['transmittance'], absorption['absorbance'] = self.AP.measure_absorption_spectrum(setting['abs_average'], setting['abs_exposure'], \
                                            absorption['reference'], filter_size=setting['filter_size'], dark_correction=True, do_plot=False)

        absorption['max'], absorption['lambda_max'] = self.AP.find_max(absorption['absorbance'], analysis_range=[300,800])

        #input('pause')

        self.mgr.pump.draw_and_dispense('flow_cell_abs', 'vial', 1*ml, velocity=1000)
        #self.mgr.pump.draw_and_dispense('air', 'vial', 1*ml, velocity=1000)


        return absorption
    

        ## start pl and abs measurement with pump
    def measure_TE(self, fname, **kwargs):

        self.send_sample_to_cell(cell = 'TE')

        self.log('Measure transient emission {timestamp}')
        #TE.measure_TE(**kwargs)

        result = self.TE.measure_TE(save_filename = fname, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **kwargs)

        self.mgr.pump.draw_and_dispense('flow_cell', 'vial', 1*ml, velocity=1000)

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
        

    def do_measurements(self, fname, measurements, dilution, **kwargs):

        if  self._check_instruments(measurements) == False:
            raise Exception()

        results, uv, PL, absorption, TE = {}, None, None, None, None
        
        config, setting = self._load_abs_PL_setting(**kwargs)

        # calculation to send to flow_cell
        total_volume = max(dilution*0.1*ml, flow_cell_volume_PL)
        total_volume = dilution*0.1*ml

        dilution_actual = total_volume/(0.1*ml)
        print('dilution : {:.2f}'.format(dilution_actual))

        self.sample_dilution(total_volume, dilution_actual)

        time.sleep(1)

        #do measurements
        if 'absorption' in measurements:
            absorption = self.measure_absorption(setting)

        #TODO dilution for PL

        if 'PL' in measurements:
            uv, PL = self.measure_PL(setting)
            
        if 'TE' in measurements:
            try:
                TE = self.measure_TE(fname, **kwargs)
                os.remove('%s_TE.pkl' %fname)
                del(TE['metadata']['sample'])
                config['TE'] = copy.deepcopy(TE['metadata'])
                del(TE['metadata'])
            except:
                print('TE measurement failed')

        results = {'uv' : uv,  'absorption' : absorption, 'PL' :  PL,  'TE' : TE, 'metadata' : config}
        self.show_abspl_results(results)

        if fname:
            self.AP.save_pkl(fname + '.pkl', results)
            self.write_result_csv(fname + '_AbsPL.csv', results)
            self.abspl_result_plots(fname, results)

        self.join_images(fname, glob.glob('%s*.png' %fname))

        self.log('Wash flow cells')
        self.mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)

        for i, cell in enumerate(measurements):
            if i == 0:
                n = 3
            else: n = 2
            self.cell_wash([flow_cells[cell]['name']], 0.5*ml, n)

        #clean the flow cell
        self.log('Wash vial')
        self.vial_wash(1.5*total_volume, 4)

        self.mgr.pump.draw_and_dispense('valve', 'waste', 0.3*ml, velocity=500)

        
        return results


if __name__ == '__main__':


    #mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.5*ml, velocity=500)
    #mgr.pump.draw_and_dispense('valve', 'flow_cell_PL', 0.5*ml, velocity=500)
    
    #mgr.pump.draw_and_dispense('ACN', 'flow_cell_PL', 1*ml, velocity=500)


    #measurment_abs_PL_TE#############
    import os 

    sample_info ={
        'name' : 'DSB',
        'concentration(uM)' : 10,
        'solvent' : 'ACN'
    }

    path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200813_test_abs_PL_TE'
    #path = 'C:/Users/hottakaz/Dropbox/PythonScript/kazu/data/20200807_test_abs_PL_TE'
    if os.path.exists(path) == False:
       os.mkdir(path)
    fname = '%s/%s_%s_in_%s_ACN-water_ref_2' %(path, sample_info['concentration(uM)'], sample_info['name'], sample_info['solvent'])
    #flist = glob.glob('%s*.png' %fname)

    Opt = Optical_measurements(TE = False, AbsPL = True)
    #Opt.join_images(fname, flist)
    # Opt.do_measurements(fname, measurements = ['absorption', 'PL', 'TE'], dilution = 1, **sample_info)
    Opt.do_measurements(fname, measurements = ['absorption'], dilution = 1, **sample_info)
    ############
    #Opt.mgr.pump.draw_and_dispense('flow_cell_abs', 'vial', 1*ml, velocity=1000)
    # Opt.mgr.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)
    # Opt.cell_wash(['flow_cell_abs'], 0.5*ml, 3)
    
    # Opt.vial_wash(0.15*ml, 2)
    #Opt.mgr.pump.draw_and_dispense('valve', 'flow_cell_abs', 3*ml)
    #Opt.mgr.pump.draw_and_dispense('valve', 'vial', 0.2*ml)
    #Opt.mgr.pump.draw_and_dispense('vial', 'waste', 1*ml)
    #Opt.mgr.pump.draw_and_dispense('valve', 'waste', 3*ml)

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






