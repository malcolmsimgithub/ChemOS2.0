import numpy as np
import matplotlib.pyplot as plt
import os, glob, time, copy, math, datetime, shutil
import pickle
import json
import csv
import random
try:
    from utils.classes import Absorption, UV, Emission
    from configs.config_AutoOpt_ocean import config

except:
    from .utils.classes import Absorption, UV, Emission
    from .configs.config_AutoOpt_ocean import config


from logger import Logger
import shutil
import sys


filedir = os.path.dirname(os.path.realpath(__file__))

DUMMY_CSV = os.path.join(filedir,"simulator_storage/1_c_1_vial5_AbsPL.csv")
DUMMY_PKL = os.path.join(filedir,"simulator_storage/1_c_1_vial5.pkl")
DUMMY_PNG = os.path.join(filedir,"simulator_storage/1_c_1_vial5.png")

today = datetime.date.today()
today_formatted = today.strftime("%Y%m%d")


config_file = '%s/configs/.config_AutoOpt_ocean.dat' %filedir

log_file = '%s/log/%s.txt' %(filedir, today_formatted)

class Optical_measurements():
    def __init__(self, TE=True, AbsPL=True, Pump=True, Evap=True, logger = False, power_control = True):

        if logger:
            self.log = Logger(stdout=True, logfile=log_file, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')
        else:
            self.log = Logger(stdout=True, logfile=None, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')

        self.config = config
    
    def do_measurements(
        self, 
        fname, 
        measurements, 
        sample_position, 
        sample_volume, 
        dilution, solvent, 
        water_ratio, 
        fill_ref = False, 
        job = None, 
        sample_info = None, 
        **kwargs):

        if job:
            self.log('\nmeasurement start for {}'.format(job['injection_name']))
        elif sample_info and sample_info['name']:
            self.log('\nmeasurement start for {}'.format(sample_info['name']))
        else :
            self.log(f'\nmeasurement start {timestamp()}')
        
        if random.uniform(0, 1) < 0.3:
            #crash the optics table
            self.log.silasocket.close()
            sys.exit()

        config = {"PL":{}, "solvent":{}}

        results = {'uv' : None,  'absorption' : None, 'PL' :  None,  'TE' : None, 'metadata' : config}

        results['metadata']['experiment_time'] = '{:%Y/%m/%d %H:%M:%S}'.format(datetime.datetime.now())

        time.sleep(2)
        if 'absorption' in measurements:
            if job:
                reference_absorption = job['average_absorbance_peak']
            else: 
                reference_absorption = None
            self.log(f'Measure absorption reference {timestamp()}')
            self.log(f'Measure absorption dark {timestamp()}')
            self.log(f'Measure absorption spectrum  {timestamp()}')

            absorption = Absorption()
            absorption.lambda_max = 123
            absorption.reference = 143
            absorption.transmittance = random.random()
            results['absorption'] = absorption._data
        
        time.sleep(2)

        if 'PL' in measurements:
            self.log(f'Measure reference uv power {timestamp()}')
            self.log(f'Adjust PL exposure {timestamp()}')
            self.log(f'Measure PL dark spectrum {timestamp()}')
            self.log(f'Measure PL spectrum {timestamp()}')
            uv, PL = UV(), Emission()

            uv.absorption = random.random()
            uv.degradation_rate= random.random()
            
            results['uv'] = uv._data
            results['PL'] = PL._data
            results['metadata']['PL']['excitation_wavelength'] = 333

        time.sleep(2)
            
        if 'TE' in measurements:
            TE = {"fitted_data": "success"}
            self.log(f'Measure transient emission {timestamp()}')
            results['TE'] = TE
        
        results['metadata']['solvent'] = "ACN"

        if job:
            results['job'] = job

        if sample_info:
            results['metadata']['sample'] = sample_info
        
        self.log('Wash cells')

        self.log(f'Wash collector {13}')

        ##create dummy pkl results file


        # with open(DUMMY_PKL, "rb") as f:
        #         dummy_pkl_data = pickle.load(f)
        # with open(f"{fname}.pkl", "wb") as f:
        #     pickle.dump(dummy_pkl_data, f)

        # ### create dummy png plots
        # shutil.copy(DUMMY_PNG, f"{fname}.png")
        # try:
        #     self.write_result_csv(f"{fname} + '_AbsPL.csv'", uv, PL, absorption)
        # except:

        #     ### create dummy csv file
        #     with open(DUMMY_CSV, "r") as f:
        #         dummy_csv_data = f.read()
        #     with open(f"{fname}_AbsPL.csv", "w") as f:
        #         f.write(dummy_csv_data)

        self.log(f'discard sample {timestamp()}')
        self.log(f'measurement done {timestamp()}')
        

        return results
    

    def auto_measurement(self, measurements = ['absorption','PL', 'TE'], solvent = 'ACN', redissolution = False, measure_blank = False):

        def blank_measurement():

            #save_dir = '%s%s' %(self.config['data_path']['result_output'], today_formatted)
            save_dir = 'output_folder/'
            if os.path.exists(save_dir) == False:
                os.mkdir(save_dir)

            print('measuring blank')
            fname = save_dir+'blank'+today_formatted+'/blank'+today_formatted
            # TODO: change
            self.measure_blank(fname, measurements = measurements, sample_position = 2, sample_volume = 0.2, \
                                            sample_info ={'name' : 'Blank', 'concentration(uM)' : 0, 'solvent' : solvent}) 
            print('blank measurement done')

        input_folder, processing_folder = os.path.abspath("input_folder/"), os.path.abspath("processing_folder/")
        completed_folder = os.path.abspath("completed_folder/")

        if measure_blank:
            blank_measurement()

        num_redisol = 0 #counter for the samples for redissolution
        waiting_list = [] #list of jobs waiting for redissolution

        self.log('sample waiting for re-dissolution %s' %num_redisol) 
        self.log(f'waiting job input {timestamp()}')

        while True:
            file_names = glob.glob(input_folder + '/*.json')
            if len(file_names) > 0:
                print('# --> file_names', file_names)
                for file_name in file_names:
                    fname = os.path.basename(file_name)
                    print(fname)
                    if fname == 'blank.pkl':
                        if 'TE' in measurements:
                            self.TE.detector_on()
                        blank_measurement()
                        os.remove(file_name)
                        continue
                    elif fname == 'shutdown.pkl':
                        self.close()
                        os.remove(file_name)
                        continue
                    else: # load parameter file
                        with open(file_name, 'r') as f:
                            file_content = json.load(f)
                            no = file_content["name"]

                        print(f"processing new job {no}")
                        save_dir = 'output_folder/'+no
                        if os.path.exists(save_dir) == False:
                            os.mkdir(save_dir)
                        fname = save_dir+"/"+file_content["name"]
                        processingname = processing_folder+"/"+no+ ".json"
                    shutil.move(file_name, processing_folder)

                    if redissolution:
                        num_redisol += 1
                        waiting_list.append([fname, file_content])

                    else:
                        self.do_measurements(fname, measurements = measurements, sample_position = file_content['vial_number'], \
                            sample_volume = file_content['sample_volume'], dilution = 1,  solvent = solvent, water_ratio = 0, sample_info = None) 
                        
                        print(processingname)
                        # remove parameter file
                        shutil.move(processingname, completed_folder+"/"+no+".json")

                        time.sleep(3)
                    
        

                input_time = time.time()
                self.log('sample waiting for re-dissolution %s' %num_redisol) 
                self.log(f'waiting job input {timestamp()}')

            if redissolution:
                if num_redisol >= self.config['redissolution']['N_parallel'] or float(time.time()-input_time) > self.config['redissolution']['waiting_timeout']:
                    self.log(f'Start redissolution process {timestamp()}')
                    self.do_redissolution([job[1]['vial_number'] for job in waiting_list])
                    # self.do_dissolution([10], solvent_volume=0.15*ml, temparature=40, duration = 300, mixing_speed=8, after_temp=25, wash_line = True, wait = True)
                    for job in waiting_list:
                        print(job[1])
                        self.do_measurements(job[1]['filename'], measurements = measurements, sample_position = job[1]['vial_number'], \
                            sample_volume = self.config['redissolution']['solvent_volume'], dilution = 1, solvent = solvent, water_ratio = 0, job = job[1]) 

                    num_redisol = 0 #reset counter
                    waiting_list = [] #reset waiting list
                    os.remove(processing_folder + job[0]) #remove job file
                    time.sleep(3)

            time.sleep(5)
        
    
    def measure_blank(self, fname, measurements, sample_position, sample_volume = 0.2, solvent = 'ACN', **kwargs):

        self.log(f'washing selector to pump line with solvent {timestamp()}')      

        # self.selector.move_to(6) #air

        self.log(f'Discard selector to pump volumes {timestamp()}')

        #self.update_vial_status(sample_position, 'used')
        time.sleep(3)

        self.do_measurements(fname, measurements = measurements, sample_position = sample_position, \
                   sample_volume = sample_volume, dilution = 1, solvent = solvent, water_ratio = 0, **kwargs) 

    

    def _load_setting(self, **kwargs):
        import copy
        config = copy.deepcopy(self.config)

        for val in config.values():
            for key in val.keys():
                if key in kwargs.keys():
                    val[key] = kwargs[key]
        
        #setting = ConfigParser(config)
        #setting.parse()

        #return config, setting
    
    def write_result_csv(self, fname, uv, PL, absorption):

        if absorption == None and PL == None:
            return

        header1, header2, header3, stats, d = [], [] ,[], [], []

        if  absorption:
            header1.extend(['abs_max', 'abs_lambda_max', 'abs_end'])
            stats.extend([absorption.max, absorption.lambda_max, absorption.end_wavelength])

            header2.extend(['Absorption', '', '', '', '', ''])
            header3.extend(['wavelength/nm', 'abs_ref', 'abs_sample', 'transmittance', 'absorbance', ''])
            d.extend([absorption.reference[0], absorption.reference[1], absorption.sample[1], \
                absorption.transmittance[1], absorption.absorbance[1], [None for _ in range(len(absorption.reference[0]))]]) 
            
        if PL:
            header1.extend(['uv_ref(W)', 'uv_absorption(W)', 'uv_absorbance_maintenance','uv_absorbance_maintenance(at_1min)', \
                                            'degradation_rate(s-1)', 'PL_max','PL_lambda_max', 'relative_QY', 'max_gain_factor(cm2 s)'])
            stats.extend([uv.reference, uv.absorption, uv.absorbance_maintenance, uv.absorbance_maintenance_at_1min,\
                                            uv.degradation_rate, PL.max, PL.lambda_max, PL.relative_QY, PL.max_gain_factor])

            header2.extend(['PL','','','','','', 'PL time trace (photons/s/nm), time/s'])
            header3.extend(['wavelength/nm', 'PL(energy/s/nm)', 'PL(photons/s/nm)', 'PL(photons/s/Hz)', 'gain_factor(cm2 s)', '', 'wavelength/nm'])
            header3.extend(PL.int_time_trace[0])
            d.extend([PL.energy[0], PL.energy[1], PL.photons[1], PL.freq_spectrum[1], PL.gain_spectrum[1], [None for _ in range(len(PL.energy[0]))], PL.energy[0]])
            d.extend(PL.time_trace) 

        with open(fname, 'w', newline = '') as f:

            writer = csv.writer(f)
            writer.writerow(header1)
            writer.writerow(stats)
            writer.writerow([])
            writer.writerow(['uv absorbance time trace'])
            writer.writerow(['time/s'] + uv.absorbance_time_trace[0].tolist())
            writer.writerow(['absorbance'] + uv.absorbance_time_trace[1].tolist())
            writer.writerow([])
            d = np.transpose(np.asarray(d))
            writer.writerow(header2)
            writer.writerow(header3)
            writer.writerows(d)
        
        print('%s was saved' %fname)
    
    #for auto measurements############################
    def load_file_content(self, file_name):
        with open(file_name, 'rb') as content:
            return pickle.load(content)
    
def timestamp():
    return time.strftime("%y-%m-%d-%H-%M", time.localtime())




        









