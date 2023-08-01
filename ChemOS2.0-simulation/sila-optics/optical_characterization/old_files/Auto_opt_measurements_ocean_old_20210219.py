from pylab.manager import Logger
import numpy as np
import matplotlib.pyplot as plt
import os, glob, time, copy, math, datetime, shutil
import pickle
import json
import csv
from pylab.instruments import PumpPSD8
from pylab.instruments import Valco_valve
from pylab.instruments import TorreyPinesSC20
from pylab.instruments import Numato_Usbrelay
from pylab.instruments import ThorlabsDC4100
from pylab.instruments import DH_mini

try:
    from transient_emission import Transient_Emission
    from absorption_and_PL_ocean import Abs_PL
    from utils.proc_image import ImageHandler
    from utils.classes import Absorption, UV, Emission
    from configs.config_parser import ConfigParser
    from power_manager.power_manager import PowerManager
    from configs.config_AutoOpt_ocean import refractive_indices, config

except:
    from .transient_emission import Transient_Emission
    from .absorption_and_PL_ocean import Abs_PL
    from .utils.proc_image import ImageHandler
    from .utils.classes import Absorption, UV, Emission
    from .configs.config_parser import ConfigParser
    from .power_manager.power_manager import PowerManager
    from .configs.config_AutoOpt_ocean import refractive_indices, config


"""
Automatic optical measurement with flow selector
"""

today = datetime.date.today()
today_formatted = today.strftime("%Y%m%d")

filedir = os.path.dirname(os.path.realpath(__file__))
config_file = '%s/configs/.config_AutoOpt_ocean.dat' %filedir

log_file = '%s/log/%s.txt' %(filedir, today_formatted)


## settings
ml = 1e-3

buffer_volume = 0.025*ml
loop_to_pump  = 0.04*ml # (under)
selector_to_pump  = 0.05*ml # (under)

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
        'to_cell_under' : 0.33*ml,
        'to_cell_over' : 0.38*ml
    } 
}

for val in flow_cells.values():
    val['cell_volume'] = round(val['to_cell_over'] - val['to_cell_under'] + 2*buffer_volume, 6)


pump_setting =   { 
            "visa": "ASRL2::INSTR",
            "addr": 0,
            "syringe_volume": 1E-3,
            "init_valve": 7,
            "ports": {
                "vial": 1,
                "ACN": 2,
                "flow_cell_PL": 3,
                "flow_cell_abs" : 4,
                "flow_cell": 5,
                "air":6,
                "waste": 7,
                "valve": 8,
            }
}


def power_on(measurements): #AbsPL and/or TE
    PM = PowerManager()
    for key, val in PM.dev_list.items():
        if val['measurement'] in measurements:
            PM.set_state(key, 'on')

def power_off(measurements): #AbsPL and/or TE
    PM = PowerManager()
    for key, val in PM.dev_list.items():
        if val['measurement'] in measurements:
            PM.set_state(key, 'off')

def lamp_on():
    DM = DH_mini('ASRL10::INSTR', verbose=False)
    DM.halogen_on()
    DM.deuterium_on()

def lamp_off():
    DM = DH_mini('ASRL10::INSTR', verbose=False)
    DM.halogen_off()
    DM.deuterium_off()

def led_on(percentage):
    led_driver = ThorlabsDC4100("ASRL9::INSTR")
    led = led_driver.led[1]
    led.set_percentage(percentage)
    led.on()

def led_off():
    led_driver = ThorlabsDC4100("ASRL9::INSTR")
    led = led_driver.led[1]
    led.off()

class Optical_measurements():

    def __init__(self, TE=True, AbsPL=True, Pump=True, Evap=True, logger = False, power_control = True):

        self.device = {'TE' : TE, 'AbsPL':AbsPL, 'Pump':Pump, 'Evap':Evap}

        if logger:
            self.log = Logger(stdout=True, logfile=log_file, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')
        else:
            self.log = Logger(stdout=True, logfile=None, pause=False, time_format='({:%Y-%m-%d %H:%M:%S})')

        # self.config = self._load_config_dat(config_file)
        self.config = config

        if self.device['Pump']:
            self.pump = PumpPSD8(**pump_setting)
            self.pump.set_velocity(1000)

            self.update_valve_status('free')

        if power_control:  #turn on devices necessary for the measurements
            power_on([key for key in self.device.keys() if self.device[key] == True])

        self.AP = Abs_PL(led_power=self.config['PL']['led_power'], correction_datafile= self.config['data_processing']['correction_datafile'], device=self.device['AbsPL'])
        self.TE = Transient_Emission(device = self.device['TE'], DB = False)
        if TE:      
            self.TE.detector_on()        
        if Evap:
            self.shaker = TorreyPinesSC20('ASRL21::INSTR') 
            self.nitrogen = Numato_Usbrelay('ASRL14::INSTR')

            self.shaker.set_temparature(self.config['redissolution']['idling_temparature'])
            self.shaker.mixer_off()

        self.IM = ImageHandler()

    def __del__(self):
        if self.device['AbsPL']:
            self.AP.led_driver.manager.close()
            self.AP.powermeter.manager.close()
            # self.AP.spectrometer.manager.close()
            self.AP.lamp.gpio.manager.close()
        if self.device['Pump']:
            self.pump.manager.close()
        if self.device['TE']:
            self.TE.fw.manager.close()
        

    def close(self):
        self.log('closing device')
        if self.TE.device:
            self.TE.detector_off()
            self.TE.laser_off()
        if self.AP.device:
            self.AP.led_off()
            self.AP.lamp_off()
        #self.log(self.AP.lamp.check_status())
        self.log(self.AP.led_365.get_state())



    #load setting#########################
    def _load_config_dat(self, config_file):
        
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
        self.log('load %s' %config_file)
        return config


    def _load_setting(self, **kwargs):
        import copy
        config = copy.deepcopy(self.config)

        for val in config.values():
            for key in val.keys():
                if key in kwargs.keys():
                    val[key] = kwargs[key]
        
        setting = ConfigParser(config)
        setting.parse()

        return config, setting


   #deal with vial status#################################
    def update_vial_status(self, vial_number, status):
        fname = '%s/available_vials.pkl' %self.config['data_path']['status']
        with open (fname, 'rb') as f:
            collector_status = pickle.load(f)
        collector_status[str(vial_number)] = status
        with open (fname, 'wb') as f:
            pickle.dump(collector_status, f)
        self.log('update vial status : vial%s %s' %(vial_number, status))


    def check_vial_status(self):
        fname = '%s/available_vials.pkl' %self.config['data_path']['status']
        with open (fname, 'rb') as f:
            print(pickle.load(f))


    #handle valve connection#################################
    def update_valve_status(self, status):
        fname = '%s/valve_status_characterization.pkl' %self.config['data_path']['status']
        with open (fname, 'wb') as f:
            pickle.dump(status, f)


    def check_valve_status(self):
        fname = '%s/valve_status_HPLCMS.pkl' %self.config['data_path']['status']
        with open (fname, 'rb') as f:
            status = pickle.load(f)
            print('HPLCMS_valve_status : %s' %status)
        return status


    def open_valves(self):
        while True:
            if self.check_valve_status() == 'free':
                self.valve = Valco_valve('visa://10.22.1.20/ASRL4::INSTR', dev_id = 0, mode = 1, position = 'A')
                self.selector = Valco_valve('visa://10.22.1.20/ASRL17::INSTR', dev_id = 0, mode = 3, position = 1)
                self.update_valve_status('busy')
                self.log('valve connection opened')
                break
            else : 
                print('Valves are busy for HPLCMS. Wait a moment\n***********************')
                time.sleep(5)


    def close_valves(self):
        # self.valve.move_to('B')
        self.selector.move_to(1)
        self.valve.close_device()
        self.selector.close_device()
        self.update_valve_status('free')
        del self.selector, self.valve
        self.log('valve connection closed')

    #overwrite pump methods###################
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

        while volume >= self.pump.syringe_volume:
            volume -= self.pump.syringe_volume
            self.draw_full(from_valve, draw_velocity)
            self.dispense_all(to_valve, dispense_velocity)

        self.draw(from_valve, volume, draw_velocity)
        self.dispense(to_valve, volume, dispense_velocity)


    def measure_line(self, line, direction='dispense', volume=1*ml, step = 0.01*ml):
        #direction : draw or dispense
        if direction == 'dispense':
            self.draw_and_dispense('air', line, 1*ml, draw_velocity =1000, dispense_velocity = 1000)
            self.draw('ACN', volume, 1000)
            action = self.pump.dispense
        else:
            action = self.pump.draw

        c = ''
        volume = 0
        while c != 'c':
            action(volume=step, valve=line)
            volume += step *1000
            c = input(f'{volume:.3f}')

        self.pump.dispense_all('waste')

    #system washing functions#################################
    def cell_wash(self, cells, volume, repeat):
        for cell in cells:
            self.log('Wash %s' %cell)
            for _ in range(repeat):
                self.draw_and_dispense('ACN', cell, volume, draw_velocity=500, dispense_velocity=1000)


    def vial_wash(self, volume, repeat, fill = False):
        self.log('Wash vial')
        for _ in range(repeat):
            self.draw_and_dispense('ACN', 'vial', volume, draw_velocity = 1000, dispense_velocity = 1000)
            self.draw_and_dispense('vial', 'waste', 1*ml, draw_velocity = 1000, dispense_velocity = 2000)
        if fill:
            self.draw_and_dispense('ACN', 'vial', volume, draw_velocity = 1000, dispense_velocity = 1000)


    def vial_empty(self):
        self.log('Send solution from the vial to waste')
        self.draw_and_dispense('vial', 'waste', 1*ml, draw_velocity = 1000, dispense_velocity = 2000)


    def collector_wash(self, volume, repeat, collector_nums, fill = False):

        self.open_valves() #open valve connection with position A and 1
        
        #cleaning the pump to selector line not to draw the solution in the line to the pump
        self.draw_and_dispense('ACN', 'valve', 0.3*ml,  draw_velocity = 1000, dispense_velocity = 1000)
        self.draw_and_dispense('air', 'valve', 0.5*ml,  draw_velocity = 1000, dispense_velocity = 1000)   

        self.valve.move_to('A')

        for collector_num in collector_nums:
            self.log(f'Wash collector {collector_num:.0f}')
            for _ in range(repeat):
                self.selector.move_to(collector_num)
                self.draw_and_dispense('ACN', 'valve', volume, draw_velocity = 1000, dispense_velocity = 1000)
                self.draw_and_dispense('valve', 'waste', 1*ml, draw_velocity = 1000, dispense_velocity = 2000)
            if fill:
                self.draw_and_dispense('ACN', 'valve', volume, draw_velocity = 1000, dispense_velocity = 1000)
            else:
                self.draw_and_dispense('valve', 'waste', 1*ml, draw_velocity = 2000, dispense_velocity = 2000)

            self.update_vial_status(collector_num, True)

        self.selector.move_to(1)
        self.close_valves()


    def collector_empty(self, collector_nums):

        self.open_valves()  #open valve connection with position A and 1

        #cleaning the pump to selector line not to draw the solution in the line to the pump
        self.draw_and_dispense('ACN', 'valve', 0.3*ml,  draw_velocity = 1000, dispense_velocity = 1000)
        self.draw_and_dispense('air', 'valve', 0.5*ml,  draw_velocity = 1000, dispense_velocity = 1000)   

        self.valve.move_to('A')

        for collector_num in collector_nums:
            self.selector.move_to(collector_num)
            self.log(f'Send solution in collector {collector_num:.0f} to waste')
            self.draw_and_dispense('valve', 'waste', 2*ml, draw_velocity = 1000, dispense_velocity = 1000)

        self.selector.move_to(1)
        self.close_valves()


    def change_solvent(self, cells):
        """  
        Replace system with different solvent.

        Parameters
        --------------------
        cells : list of str
            cells to be used. (current options : 'TE', 'PL', 'absorption')
        """
        cell_names = [flow_cells[key]['name'] for key in flow_cells.keys() if key in cells]

        #first discard all previous solutions to waste 
        for i, cell in enumerate(cell_names):
            self.log('Discard solution in %s' %cell)
            self.draw_and_dispense(cell, 'waste', 2*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
        #then replace the syringe dead volume to new solvent
        self.log('replacing the syringe dead volume to the new solvent')
        self.draw_and_dispense('ACN', 'waste', 2*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
        #then fill solutions to cells
        for i, cell in enumerate(cell_names): 
            self.log('sending the new solvent to %s' %cell)  
            self.draw_and_dispense('ACN', cell, 3*ml, draw_velocity = 1000, dispense_velocity = 1000)

        self.vial_wash(0.5*ml, 3)

        self.collector_wash(0.5*ml, 3, [2,3,4,5,7,8,9,10])


    #Evaporation and re-dissolution######################
    def do_evaporation(self, collector_nums, evap_temp, duration, after_temp, nitrogen = True, wait = True):
        """  
        Do evaporation using TorryPains SC20 with optional nitrogen down stream.

        Parameters
        --------------------
        collector_nums : list of int
            collector numbers to be evaporated.
        evap_temp : float
            temparature setting for evaporation
        duration : int
            duration for evaporation in seconds
        after_temp : float
            temparature setting for after evaporation
        nitrogen : bool
            turn on nitrogen stream during evaporation if True
        wait : bool
            wait until temparture reachs to after_temp if True
        """
        if self.device['Evap'] == False:
            raise Exception('<do_evaporation> : Evaporator not connected')

        # self.draw_and_dispense('ACN', 'valve', 0.3*ml, 500, 500)

        self.open_valves() #initialze valve at A and 1

        # self.log('washing selector to pump line with solvent {timestamp}')
        # self.draw_and_dispense('ACN', 'valve', 0.8*ml,  draw_velocity = 1000, dispense_velocity = 1000)  
        self.log('emptying selector to pump line with the air {timestamp}')
        self.draw_and_dispense('air', 'valve', 1*ml,  draw_velocity = 1000, dispense_velocity = 1000)  

        # self.selector.move_to(6) #air

        # self.log('Discard selector to pump volumes')
        # self.draw_and_dispense('valve', 'waste', 0.5*ml,  draw_velocity = 1000, dispense_velocity = 1000)

        #make sure all solutions in the tubing go to the collection vials

        for i in collector_nums:
            self.log('Sending solutions to the collector vials %s' %i)
            self.selector.move_to(i)
            self.draw_and_dispense('air', 'valve', 2*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
        self.close_valves() #close valves at B and 1

        # input("pause")

        self.shaker.mixer_off()
        self.shaker.set_temparature(evap_temp, wait = False)

        if nitrogen:
            self.nitrogen.relay_on(0)
            self.log('nitrogen_valve: Open')
        
        self.log('Waiting {:.1f} min for evaporation'.format(duration/60))
        self.log('{timestamp}')
        time.sleep(duration*0.8)

        # input('check the tubing')

        self.open_valves()

        self.draw_and_dispense('air', 'valve', 1*ml,  draw_velocity = 1000, dispense_velocity = 1000)  

        for i in collector_nums:
            self.log('Sending remaining solutions to the collector vials %s' %i)
            self.selector.move_to(i)
            self.draw_and_dispense('air', 'valve', 1*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
        self.close_valves() #close valves at B and 1        
        
        time.sleep(duration*0.2)

        if nitrogen:
            self.nitrogen.relay_off(0)
            self.log('nitrogen_valve: Close')
        self.shaker.set_temparature(after_temp, wait = wait, tolerance = 2)

    
    def do_dissolution(self, collector_nums, solvent_volume, temparature, duration, mixing_speed, after_temp, wash_line = False, wait = True):

        if self.device['Evap'] == False:
            raise Exception('<do_dissolution> : Evaporator not connected')

        self.shaker.set_temparature(temparature, wait = True)

        # self.draw_and_dispense('ACN', 'valve',  0.3*ml, 500, 500)

        self.open_valves()  #initialze valve at A and 1

        if wash_line:
            self.log('washing selector to pump line with solvent {timestamp}')
            # self.draw_and_dispense('air', 'valve', 0.6*ml,  draw_velocity = 1000, dispense_velocity = 1000)  
            self.draw_and_dispense('ACN', 'valve', 0.5*ml,  draw_velocity = 1000, dispense_velocity = 1000)  
            self.draw_and_dispense('air', 'valve', 1*ml,  draw_velocity = 1000, dispense_velocity = 1000) 
            # self.selector.move_to(6) #air

            # self.log('Discard selector to pump volumes')
            # self.draw_and_dispense('valve', 'waste', 0.3*ml,  draw_velocity = 1000, dispense_velocity = 1000)

        for i in collector_nums:
            self.log('Sending %s ml of solvent to collector vial %s' %(solvent_volume*1000, i)) 
            self.selector.move_to(i)
            # self.draw_and_dispense('ACN', 'waste', 0.2*ml, draw_velocity = 1000, dispense_velocity = 1000)
            self.draw_and_dispense('ACN', 'valve', solvent_volume, draw_velocity = 1000, dispense_velocity = 1000)
            self.draw_and_dispense('air', 'valve', 1*ml, draw_velocity = 1000, dispense_velocity = 1000)
        
        self.close_valves() #close valves at B and 1

        self.shaker.mixer_on(mixing_speed)

        self.log('Waiting {:.1f} min for dissolution'.format(duration/60))
        self.log('{timestamp}')
        time.sleep(duration)
        
        self.shaker.set_temparature(after_temp, wait = wait, tolerance = 2)
        self.shaker.mixer_off()


    #transfer solutions#################################
    def fill_reference(self, water_ratio, measurements):

        num = len([measurement for measurement in measurements if measurement == 'absorption' or measurement == 'PL'])

        if water_ratio == 0 or num == 0:
            pass
        else:
            self.open_valves()
            self.log('mixing reference solution')
            self.draw_and_dispense('ACN', 'vial', 0.5 * num * (1-water_ratio)*ml, draw_velocity = 1000, dispense_velocity = 1000)
            #draw water
            self.selector.move_to(1)
            self.valve.move_to('B')
            self.draw_and_dispense('valve', 'waste', 0.3*ml, draw_velocity=1000, dispense_velocity=1000) #discard tube volume
            self.draw_and_dispense('valve', 'vial', 0.5 * num * water_ratio*ml, draw_velocity = 1000, dispense_velocity = 300)
            self.valve.move_to('A')

            for _ in range(2): #mixing
                self.draw_and_dispense('vial', 'vial', 1*ml, draw_velocity = 300, dispense_velocity = 300)
            
            if 'absorption' in measurements:
                self.log('Sending reference solution to absorption cell')
                self.draw_and_dispense('vial', 'flow_cell_abs', 0.5*ml, draw_velocity = 300, dispense_velocity = 300)
            if 'PL' in measurements:
                self.log('Sending reference solution to PL cell')
                self.draw_and_dispense('vial', 'flow_cell_PL', 0.5*ml, draw_velocity = 300, dispense_velocity = 300)
            self.close_valves() #close valves at B and 1


    def collect_sample_from_valve(self, sample_volume, collector_num, setting):

        self.open_valves() #initialze valve at A and 1

        self.selector.move_to(1)

        self.log('washing selector to pump line with solvent {timestamp}')
        self.draw_and_dispense('ACN', 'valve', 0.3*ml, setting.dilution_draw_velocity, setting.dilution_dispense_velocity)        

        # self.selector.move_to(6) #air

        self.log('Discard selector to pump volumes')
        # self.draw_and_dispense('valve', 'waste', 0.3*ml, setting.dilution_draw_velocity, setting.dilution_dispense_velocity)
        self.draw_and_dispense('air', 'valve', 0.5*ml, setting.dilution_draw_velocity, setting.dilution_dispense_velocity)

        self.log(f'Switch selector to position {collector_num:.0f}')
        self.selector.move_to(collector_num)

        self.log(f'Draw {sample_volume*1000:.3f} ml of sample')
        self.draw_and_dispense('valve', 'vial', 1*ml, setting.dilution_draw_velocity, setting.dilution_dispense_velocity)
        
        #######
        # self.draw('valve', 1*ml, setting['dilution_draw_velocity'])
        # input('pause')
        # self.dispense_all('waste', setting['dilution_dispense_velocity'])
        #########
        
        #self.draw_and_dispense('valve', 'vial', 0.5*ml, setting.dilution_draw_velocity, setting.dilution_dispense_velocity)

        self.log(f'Switch selector to position 1')  

        self.close_valves() #close valves at B and 1


    def sample_dilution(self, dilution, setting):

        if dilution <= 1:
            raise Exception('dilution should be larger than 1')

        self.log(f'Draw {self.sample_volume*(dilution - 1)*1000:.3f} ml of solvent and mix in vial')
        self.draw_and_dispense('ACN', 'vial', self.sample_volume*(dilution - 1), setting.dilution_draw_velocity, setting.dilution_dispense_velocity)
        for _ in range(2):
            self.draw_and_dispense('vial', 'vial', max(self.sample_volume * dilution *4, 0.5*ml), setting.dilution_draw_velocity, setting.dilution_dispense_velocity) # mix 3 times

        self.sample_volume *= dilution
        self.log(f'Current sample volume : {self.sample_volume/ml:.3f} ml')


    def send_sample_to_cell(self, cell, draw_velocity = 500, dispense_velocity = 1000):

        flow_cell = flow_cells[cell]
        volume = flow_cell['cell_volume']

        #make flow cell empty
        self.log('Emplying the flow cell {timestamp}')
        self.pump.draw_and_dispense('air', flow_cell['name'], 0.5*ml, velocity=1000) #better not to draw solvent to the pump
        self.pump.draw_and_dispense('air', flow_cell['name'], 1*ml, velocity=2000) #to prevent unintentioanl dilution
        #self.pump.draw_and_dispense(flow_cell['name'], 'waste', 1*ml, velocity=2000) #to prevent unintentioanl dilution

        self.log(f'Send {volume/ml:.3f} ml of solution to flow cell')
        draw_volume = min(1*ml, volume * 3)
        # draw_volume = 0.15*ml
        self.draw('vial', draw_volume , draw_velocity)
        # self.draw('ACN', draw_volume , draw_velocity)

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


    #process measurement results#################################
    def save_data(self, fname, data):
        np.savetxt(fname, data, delimiter=',', header='wl (nm), abs (%), pl (and abs/pl pair repeats)')


    #TODO time trace
    def write_result_csv(self, fname, uv, PL, absorption):

        if absorption == None and PL == None:
            return

        header1, header2, stats, d = [], [] ,[], []

        if  absorption:
            header1.extend(['abs_max', 'abs_lambda_max', 'abs_end'])
            stats.extend([absorption.max, absorption.lambda_max, absorption.end_wavelength])

            header2.extend(['wavelength/nm', 'abs_ref', 'abs_sample', 'transmittance', 'absorbance'])
            d.extend([absorption.reference[0], absorption.reference[1], absorption.sample[1], absorption.transmittance[1], absorption.absorbance[1]]) 
            
        if PL:
            header1.extend(['uv_ref(W)', 'uv_absorption(W)', 'uv_absorbance_maintenance','uv_absorbance_maintenance(at_1min)', 'degradation_rate(s-1)', 'PL_max','PL_lambda_max', 'relative_QY', 'max_gain_factor(cm2 s)'])
            stats.extend([uv.reference, uv.absorption, uv.absorbance_maintenance, uv.absorbance_maintenance_at_1min, uv.degradation_rate, PL.max, PL.lambda_max, PL.relative_QY, PL.max_gain_factor])

            header2.extend(['wavelength/nm', 'PL(energy/s/nm)', 'PL(photons/s/nm)', 'PL(photons/s/Hz)', 'gain_factor(cm2 s)'])
            d.extend([PL.energy[0], PL.energy[1], PL.photons[1], PL.freq_spectrum[1], PL.gain_spectrum[1]]) 

        with open(fname, 'w', newline = '') as f:

            writer = csv.writer(f)
            writer.writerow(header1)
            writer.writerow(stats)
            d = np.transpose(np.asarray(d))
            writer.writerow(header2)
            writer.writerows(d)
        
        print('%s was saved' %fname)



    def abspl_result_plots(self, fname, uv, PL, absorption):
        #absortpion spectrum
        if absorption:
            text = 'abs_max : {:.4f} ({:.1f} nm)'.format(absorption.max, absorption.lambda_max)
            if absorption.end_wavelength:
                text = text + '\nabs_end : {:.1f} nm'.format(absorption.end_wavelength)
            self.AP.result_plot(absorption.absorbance, 'Absorbance', text=text, xrange = [200,800], show_plot=False, save_filename='%s_absorption_spactrum.png' %fname)
        #PL spectrum
        if PL:
            text = 'PL_max : {:.4f} ({:.1f} nm)\nabsorbed_power : {:.3f} mW\nrelative_QY : {:.3f}'.format\
                           (PL.max, PL.lambda_max, uv.absorption*1000, PL.relative_QY)
            if uv.absorbance_maintenance:
                text = text + '\nuv_absorbance_maintenance(%) : {:.2f}\nuv_absorbance_maintenance_1min(%) : {:.2f}\ndegradation_rate(s-1) : {:.5f}'.format\
                            (uv.absorbance_maintenance*100, uv.absorbance_maintenance_at_1min*100, uv.degradation_rate)

            self.AP.result_plot(PL.energy, 'Intensity/ (energy/s/nm)', text = text, xrange = [300,800], show_plot=False, save_filename='%s_PL_spactrum.png' %fname)
            self.AP.result_plot(PL.gain_spectrum, 'gain factor/ (cm2 s)', text = None, xrange = [300,800], show_plot=False, save_filename='%s_PL_gain_spactrum.png' %fname)
            self.AP.plot_time_trace(uv._data, PL._data, show_plot=False, save_filename='%s_time_trace.png' %fname)
        #abs/PL spectrum
        if absorption and PL:
            self.AP.Abs_PL_plot(absorption.absorbance, PL.energy, text = None, xrange = [300,800], show_plot=False, save_filename='%s_Abs_PL_spactrum.png' %fname)
            if PL.max_gain_factor:
                text = 'max_gain_factor(cm2 s) :\n {:.3e} ({:.1f} nm)'.format(PL.max_gain_factor, PL.max_gain_wavelength)
            else:
                text = None
            self.AP.Abs_PL_plot(absorption.absorbance, PL.gain_spectrum, text = text, xrange = [300,800],\
                                            ylabel = 'Nomalized Abs. and gain_factor.', show_plot=False, save_filename='%s_Abs_PL_gain_spectrum.png' %fname)



    def PL_uv_analysis(self, PL_uv_results, uv_reference, PL_exposure, calc_range, cut_off = 0.95):

        l_index, u_index = self.AP._to_index_range(self.AP.wl, calc_range)
        uv = UV() 
        PL = Emission()

        uv.reference = uv_reference
        uv.absorbance_time_trace = np.asarray([PL_uv_results['uv_time_trace'][0]-PL_uv_results['uv_time_trace'][0, 0], -np.log10(PL_uv_results['uv_time_trace'][1]/uv.reference)])
        PL.time_trace = PL_uv_results['PL']['time_trace']
        PL.int_time_trace = np.asarray([PL_uv_results['PL']['time'], [np.sum(trace[l_index:u_index+1]) for trace in PL.time_trace]])
        # PL.int_time_trace = np.asarray([PL_uv_results['PL']['time'], [np.mean(trace[max_index-3:max_index+3]) for trace in PL.time_trace]])
        
        #find the time range for calculation based on the PL intensity cutoff
        PL_indice = np.where(PL.int_time_trace[1]/PL.int_time_trace[1][0] > cut_off)[0]
        cut_off_time = PL.int_time_trace[0][PL_indice[-1]]
        uv_indice = np.where(uv.absorbance_time_trace[0] < cut_off_time)[0]
        if len(uv_indice) == 0:
            uv_indice = [0]

        PL_average = np.mean(PL.time_trace[PL_indice], axis = 0)
        uv.absorption = uv.reference - np.mean(PL_uv_results['uv_time_trace'][1][uv_indice])

        PL.photons = np.asarray([PL_uv_results['PL']['wavelength'], PL_average/PL_exposure])   #the value is calibrated against photons/nm
        PL.energy = np.asarray([PL_uv_results['PL']['wavelength'], PL.photons[1]/PL.photons[0]/1E-6])
        PL.freq_spectrum = self.AP.calc_freq_spectrum(PL.photons, calc_range= calc_range)
        PL.gain_spectrum = self.AP.calc_gain_spectrum(PL.photons, calc_range= calc_range)
    
        area = np.sum(PL.photons[1][l_index:u_index+1] * (PL.photons[0][l_index+1:u_index+2] - PL.photons[0][l_index-1:u_index])/2)

        PL.relative_QY = self.qy_ref['QY'] / self.qy_ref['value'] * (area/uv.absorption)  * (refractive_indices[self.solvent]/refractive_indices[self.qy_ref['solvent']]) **2 
        # PL.relative_QY = self.qy_ref['QY'] * (np.sum(PL.photons[1][l_index:u_index+1])/uv.absorption) / self.qy_ref['value'] * (refractive_indices[self.solvent]/refractive_indices[self.qy_ref['solvent']]) **2 
                                                            
        PL.max, PL.lambda_max = self.AP.find_max(PL.energy, analysis_range=calc_range)
        # max_index = np.argmax(PL.energy[1][l_index:u_index+1])+l_index

        print('PL_max :  {:.4f} ({:.1f} nm)'.format(PL.max, PL.lambda_max))
        print('uv_ref : {:.5f} mW'.format(uv.reference * 1000))
        print('uv_start : {:.5f} mW'.format(PL_uv_results['uv_time_trace'][1, 0] * 1000))
        print('uv_end : {:.5f} mW'.format(PL_uv_results['uv_time_trace'][1, -1] * 1000))
        print('relative_QY : {:.3f}'.format(PL.relative_QY))

        if uv.absorption/uv.reference > 0.001:
            uv.absorbance_maintenance = min(1, math.log10(PL_uv_results['uv_time_trace'][1, -1]/uv.reference) /math.log10(PL_uv_results['uv_time_trace'][1, 0]/uv.reference))
            uv.degradation_rate = - math.log(uv.absorbance_maintenance)/PL_uv_results['duration']
            uv.absorbance_maintenance_at_1min = uv.absorbance_maintenance * math.exp(-(60-PL_uv_results['duration'])* uv.degradation_rate)

            print('uv_absorbance_maintenance : {:.3f}%'.format(uv.absorbance_maintenance*100))
            print('uv_degradation_rate: {:.3f}'.format(uv.degradation_rate))
            print('uv_absorbance_maintenance(at_1min) : {:.3f}%'.format(uv.absorbance_maintenance_at_1min*100))

        if PL.int_time_trace[1, -1] > 0 and PL.int_time_trace[1, 0] > 0:
            PL.maintenance = min(1, PL.int_time_trace[1, -1]/PL.int_time_trace[1, 0])
            PL.degradation_rate = - math.log(PL.maintenance)/PL_uv_results['duration']
            PL.maintenance_at_1min = PL.maintenance * math.exp(-(60-PL_uv_results['duration'])* PL.degradation_rate)

            print('PL_maintenance : {:.3f}%'.format(PL.maintenance*100))
            print('PL_degradation_rate: {:.3f}'.format(PL.degradation_rate))
            print('PL_maintenance(at_1min) : {:.3f}%'.format(PL.maintenance_at_1min*100))

        return uv, PL


    def _find_maximum_gain(self, gain_spectrum , absorption_end_index, calc_range = [300,800]):

        l_index, u_index = self.AP._to_index_range(gain_spectrum[0], calc_range)
        l_index = max(l_index, absorption_end_index)

        spectrum = gain_spectrum[:, l_index:u_index+1]

        max_gain = np.max(spectrum[1])
        max_gain_wavelength = spectrum[0][np.argmax(spectrum[1])]

        print('maxinum_spectral_gain_factor : {:.3f} 1E-24 cm2 s({:.1f} nm)'.format(max_gain*1E24, max_gain_wavelength))

        return max_gain, max_gain_wavelength


    def join_images(self, fname, filelist):
        if len(filelist) > 1:
            self.IM.tile_img(h_num = 3, save_filename = '%s.png' %fname, mergin = [0,0,0,0], del_files=True, filelist=filelist)


    ###manual measurements#########################################
    def measure_PL_manual(self, fname, sample_info = None, **kwargs):

        config, setting = self._load_setting(**kwargs)

        results = {'uv' : None,  'absorption' : None, 'PL' :  None,  'TE' : None, 'metadata' : config}

        input('Please inject reference to the cell')

        #measure reference
        self.log('Measure reference uv power {timestamp}')
        uv_reference = self.AP.measure_uv_power(counts = setting.PL_uv_average, led_power = setting.PL_led_power)


        input('Please inject the sample to the cell')

        self.log('Adjust PL exposure {timestamp}')
        PL_exposure = self.AP.adjust_PL_exposure(setting.PL_initial_exposure, setting.PL_max_exposure, setting.PL_target_intensity,\
                                                                        setting.PL_led_power, filter_size=setting.filter_size, average = 5)
        if PL_exposure * setting.PL_average < setting.PL_min_measurement_time:
            PL_average = int(setting.PL_min_measurement_time/PL_exposure)
        else:
            PL_average = setting.PL_average
        self.log('Measure PL dark spectrum {timestamp}')
        self.AP.measure_dark_spectrum(20, PL_exposure, filter_size=setting.filter_size,  do_plot=False) #dark for PL
        self.log('Measure PL spectrum {timestamp}')
        res = self.AP.measure_PL_uv(PL_average,  PL_exposure,  led_power = setting.PL_led_power,\
                                    filter_size = setting.filter_size, dark_correction = True, spectral_correction = True, do_plot = False)

        #calculations
        self.qy_ref = setting.quantum_yeild_reference
        uv, PL = self.PL_uv_analysis(res, uv_reference, PL_exposure, setting.PL_calc_range, cut_off=setting.PL_cutoff)
        PL.exposure = PL_exposure

        results['uv'] = uv._data
        results['PL'] = PL._data

        if sample_info:
            results['metadata']['sample'] = sample_info

        if fname:
            self.AP.save_pkl(fname + '.pkl', results)
            self.write_result_csv(fname + '_AbsPL.csv', uv, PL, None)
            self.abspl_result_plots(fname, uv, PL, None)

            self.join_images(fname, glob.glob('%s*.png' %fname))


        input('Measurement done. Please wash the cell')

        return uv, PL


    #measurement fuctions#################################
    def measure_PL(self, setting):

        results = {}

        #dilute the sample
        dilution = max(1, setting.PL_minimum_volume * ml/self.sample_volume)

        print('dilution for PL measurement: {:.2f}'.format(dilution))
        if dilution > 1:
            self.sample_dilution(dilution, setting)

        #measure reference
        self.log('Measure reference uv power {timestamp}')
        uv_reference = self.AP.measure_uv_power(counts = setting.PL_uv_average, led_power = setting.PL_led_power)

        self.send_sample_to_cell(cell = 'PL', draw_velocity=setting.PL_draw_velocity, dispense_velocity=setting.PL_dispense_velocity)
        # input("pause")
        # need sometime until the solution reaches equiliblium (bubble??)
        print('waiting %s s for equibliration' %setting.PL_equibliration_time)
        time.sleep(setting.PL_equibliration_time) 


        self.log('Adjust PL exposure {timestamp}')
        PL_exposure = self.AP.adjust_PL_exposure(setting.PL_initial_exposure, setting.PL_max_exposure, setting.PL_target_intensity,\
                                                                        setting.PL_led_power, filter_size=setting.filter_size, average = 5)
        if PL_exposure * setting.PL_average < setting.PL_min_measurement_time:
            PL_average = int(setting.PL_min_measurement_time/PL_exposure)
        else:
            PL_average = setting.PL_average
        self.log('Measure PL dark spectrum {timestamp}')
        self.AP.measure_dark_spectrum(20, PL_exposure, filter_size=setting.filter_size,  do_plot=False) #dark for PL
        self.log('Measure PL spectrum {timestamp}')
        results = self.AP.measure_PL_uv(PL_average,  PL_exposure,  led_power = setting.PL_led_power,\
                                    filter_size = setting.filter_size, dark_correction = True, spectral_correction = True, do_plot = False)
        # results = self.AP.measure_PL_uv(setting.PL_average,  PL_exposure,  led_power = setting.PL_led_power,\
        #                             filter_size = setting.filter_size, dark_correction = True, spectral_correction = True, do_plot = False)

        #calculations
        self.qy_ref = setting.quantum_yeild_reference
        uv, PL = self.PL_uv_analysis(results, uv_reference, PL_exposure, setting.PL_calc_range, cut_off=setting.PL_cutoff)
        PL.exposure = PL_exposure

        # input('pause')
        self.collect_sample_from_cell(cell = 'PL', draw_velocity=300, dispense_velocity=300)

        return uv, PL


    def measure_absorption(self, setting, reference_absorption = None):

        Ab = Absorption()

        #dilute the sample
        if reference_absorption:
            dilution = max(1, float(reference_absorption)/setting.abs_maximum_absorption/1000, setting.abs_minimum_volume* ml/self.sample_volume)
        else:
            dilution = max(1, setting.abs_minimum_volume* ml/self.sample_volume)

        print('dilution for absorption measurement: {:.2f}'.format(dilution))
        if dilution > 1:
            self.sample_dilution(dilution, setting)

        #measure reference
        self.log('Measure absorption reference {timestamp}')
        # self.AP.measure_dark_spectrum(setting.abs_dark_average, setting.abs_exposure, filter_size=setting.filter_size *3,  do_plot=False)
        self.AP.measure_dark_spectrum(setting.abs_dark_average, setting.abs_exposure, filter_size=setting.filter_size *3, do_plot=False)
        Ab.reference = self.AP.measure_transmission_spectrum(setting.abs_average, setting.abs_exposure, filter_size = setting.filter_size, dark_correction=True, do_plot=False)

        self.send_sample_to_cell(cell = 'absorption', draw_velocity=setting.abs_draw_velocity, dispense_velocity=setting.abs_dispense_velocity)
        # self.cell_wash([flow_cells['absorption']['name']], 0.5*ml, 1)

        #need sometime until the solution reaches equiliblium (bubble??)
        print('waiting %s s for equibliration' %setting.abs_equibliration_time)
        time.sleep(setting.abs_equibliration_time) 

        #measure absorption spectrum
        self.log('Measure absorption dark {timestamp}')
        self.AP.measure_dark_spectrum(setting.abs_dark_average, setting.abs_exposure, filter_size=setting.filter_size *3, do_plot=False)
        self.log('Measure absorption spectrum  {timestamp}')
        Ab.sample, Ab.transmittance, Ab.absorbance = self.AP.measure_absorption_spectrum(setting.abs_average, setting.abs_exposure, \
                                            Ab.reference, filter_size=setting.filter_size, dark_correction=True, do_plot=setting.abs_do_plot)

        Ab.max, Ab.lambda_max = self.AP.find_max(Ab.absorbance, analysis_range=setting.abs_calc_range) 
        Ab.end_index, Ab.end_wavelength = self.AP.find_abs_end(Ab.absorbance, setting.absorption_threshold, analysis_range=setting.abs_calc_range)     
        print('abs_max :  {:.4f} ({:.1f} nm)'.format(Ab.max, Ab.lambda_max))
        if Ab.end_wavelength:
            print('abs_end :  {:.1f} nm'.format(Ab.end_wavelength))

        #input('pause')
        self.collect_sample_from_cell(cell = 'absorption', draw_velocity=300, dispense_velocity=300)

        return Ab
    

        ## start pl and abs measurement with pump
    def measure_TE(self, fname, setting):

        #dilute the sample
        dilution = max(1, setting.TE_minimum_volume* ml/self.sample_volume)

        print('dilution for TE measurement: {:.2f}'.format(dilution))
        if dilution > 1:
            self.sample_dilution(dilution, setting)

        self.send_sample_to_cell(cell = 'TE', draw_velocity=setting.TE_draw_velocity, dispense_velocity=setting.TE_dispense_velocity)

        self.log('Measure transient emission {timestamp}')

        result = self.TE.measure_TE(save_filename = fname, do_plot = False, **setting.TE.to_dict())

        self.collect_sample_from_cell(cell = 'TE', draw_velocity=setting.TE_draw_velocity, dispense_velocity=setting.TE_dispense_velocity)
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


    def do_redissolution(self, collector_nums, **kwargs):

        _, setting = self._load_setting(**kwargs)

        self.do_evaporation(collector_nums, setting.evap_temp, setting.evap_time, setting.dissolution_temp, nitrogen=True, wait = True)

        # input("pause")

        self.do_dissolution(collector_nums, setting.dissolution_volume*ml, setting.dissolution_temp, setting.dissolution_time, \
                                                                    setting.shaking_speed, after_temp=setting.dissolution_idling_temp, wash_line=True, wait=True)

        # input("pause")


    #TODO keep sample after the measurements
    def do_measurements(self, fname, measurements, sample_position, sample_volume, dilution, solvent, water_ratio, job = None, sample_info = None, **kwargs):

        if solvent not in refractive_indices.keys():
            raise Exception('refractive index of <%s> is not specified in the config file')
        else:
            self.solvent = solvent

        if  self._check_instruments(measurements) == False:
            raise Exception()

        if job:
            self.log('\nmeasurement start for {}'.format(job['injection_name']))
        elif sample_info and sample_info['name']:
            self.log('\nmeasurement start for {}'.format(sample_info['name']))
        else :
            self.log('\nmeasurement start {timestamp}')

        config, setting = self._load_setting(**kwargs)

        self.sample_volume = sample_volume

        uv, PL, absorption, TE = None, None, None, None
        results = {'uv' : None,  'absorption' : None, 'PL' :  None,  'TE' : None, 'metadata' : config}
        
        results['metadata']['experiment_time'] = '{:%Y/%m/%d %H:%M:%S}'.format(datetime.datetime.now())
        # print(results['metadata']['experiment_time'])

        if sample_position:

            self.fill_reference(water_ratio, measurements) #TODO consider dilution on water ratio
            self.collect_sample_from_valve(self.sample_volume, sample_position, setting)

            # self.draw_and_dispense('ACN', 'vial', 0.5*ml, 500, 500)

        print('initial dilution : {:.2f}'.format(dilution))
        if dilution > 1:
            self.sample_dilution(dilution, setting)

        time.sleep(1)

        #do measurements
        if 'absorption' in measurements:
            if job:
                reference_absorption = job['average_absorbance_peak']
            else: 
                reference_absorption = None
            absorption = self.measure_absorption(setting, reference_absorption=reference_absorption)
            results['absorption'] = absorption._data

        #TODO dilution for PL

        if 'PL' in measurements:
            uv, PL = self.measure_PL(setting)
            results['uv'] = uv._data
            results['PL'] = PL._data
            
        if 'TE' in measurements:
            TE = self.measure_TE(fname, setting)
            time.sleep(1)
            if TE:
                os.remove('%s_TE.pkl' %fname)
            results['TE'] = TE
            
        if absorption is not None and PL is not None:
            if absorption.end_index:
                PL.max_gain_factor, PL.max_gain_wavelength = self._find_maximum_gain(PL.gain_spectrum, absorption.end_index, calc_range = setting.PL_calc_range)

        results['metadata']['solvent'] = self.solvent

        if job:
            results['job'] = job

        if sample_info:
            results['metadata']['sample'] = sample_info

        if fname:
            self.AP.save_pkl(fname + '.pkl', results)
            self.write_result_csv(fname + '_AbsPL.csv', uv, PL, absorption)
            self.abspl_result_plots(fname, uv, PL, absorption)
            self.join_images(fname, glob.glob('%s*.png' %fname))
        self.log('discard sample {timestamp}')
        self.pump.draw_and_dispense('vial', 'waste', 1*ml, velocity=2000)

        #wash flow cells and vial
        for i, cell in enumerate(measurements):
            if i == 0:
                n = 3
            else: n = 2
            self.cell_wash([flow_cells[cell]['name']], 0.5*ml, n)

        self.vial_wash(1.5*self.sample_volume, 3)

        if sample_position == 9 or sample_position == 10:
            self.collector_wash(0.5*ml, 3, [sample_position])
        elif sample_position:
            self.collector_wash(0.3*ml, 3, [sample_position])
        # self.collector_wash(2*self.sample_volume, 3, [sample_position])

        self.log('measurement done {timestamp}')

        return results


    #for auto measurements############################
    def load_file_content(self, file_name):
        with open(file_name, 'rb') as content:
            return pickle.load(content)


    def experiment_params(self, file_content):

        save_dir = 'C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/sample_measured/%s' %today_formatted

        if os.path.exists(save_dir) == False:
            os.mkdir(save_dir)
        file_content['filename'] = save_dir + '/{}_{}_{:.3f}min_v{}'.format(file_content['injection_name'], file_content['target_name'], file_content['retention_time'], file_content['vial_number'])

        return file_content
        

    def _check_HPLCMS_status(self):
        return os.path.exists('%s/HPLCMS_waiting.pkl' %self.config['data_path']['status'])


    def auto_measurement(self, measurements = ['absorption','PL', 'TE'], solvent = 'ACN', redissolution = False, measure_blank = False):

        def blank_measurement():

            today = datetime.date.today()
            today_formatted = today.strftime("%Y%m%d")
            save_dir = 'C:/Users/MatterLab/Dropbox/PythonScript/HPLCMS_characterization/sample_measured/%s' %today_formatted
            if os.path.exists(save_dir) == False:
                os.mkdir(save_dir)

            print('measuring blank')
            fname = '%s/blank' %save_dir
            self.measure_blank(fname, measurements = measurements, sample_position = 2, sample_volume = 0.2*ml, \
                                            sample_info ={'name' : 'Blank', 'concentration(uM)' : 0, 'solvent' : solvent}) 
            print('blank measurement done')

        input_folder, processing_folder = self.config['data_path']['job_input'], self.config['data_path']['job_processing']

        if measure_blank:
            blank_measurement()

        num_redisol = 0 #counter for the samples for redissolution
        waiting_list = [] #list of jobs waiting for redissolution

        self.log('sample waiting for re-dissolution %s' %num_redisol) 
        self.log('waiting job input {timestamp}')

        while True:
            file_names = glob.glob(input_folder + '/*pkl')
            if len(file_names) > 0:
                print('# --> file_names', file_names)
                for file_name in file_names:
                    fname = os.path.basename(file_name)
                    print(fname)
                    if fname == 'blank.pkl':
                        if 'TE' in measurements:
                            self.TE.detector_on()
                        blank_measurement()
                    elif fname == 'shutdown.pkl':
                        self.close()
                    else: # load parameter file
                        file_content = self.load_file_content(file_name)
                        file_content = self.experiment_params(file_content)
                        print(file_content)

                    shutil.move(file_name, processing_folder)

                    if redissolution:
                        num_redisol += 1
                        waiting_list.append([fname, file_content])

                    else:
                        self.do_measurements(file_content['filename'], measurements = measurements, sample_position = file_content['vial_number'], \
                            sample_volume = file_content['sample_volume'] *ml, dilution = 1,  solvent = solvent, water_ratio = 0, job = file_content) 
        
                        # remove parameter file
                        os.remove(processing_folder + fname) #remove job file
                        time.sleep(3)

                input_time = time.time()
                self.log('sample waiting for re-dissolution %s' %num_redisol) 
                self.log('waiting job input {timestamp}')

            if redissolution:
                if num_redisol >= self.config['redissolution']['N_parallel'] or self._check_HPLCMS_status() or float(time.time()-input_time) > self.config['redissolution']['waiting_timeout']:
                    self.log('Start redissolution process {timestamp}')
                    self.do_redissolution([job[1]['vial_number'] for job in waiting_list])
                    # self.do_dissolution([10], solvent_volume=0.15*ml, temparature=40, duration = 300, mixing_speed=8, after_temp=25, wash_line = True, wait = True)
                    for job in waiting_list:
                        print(job[1])
                        self.do_measurements(job[1]['filename'], measurements = measurements, sample_position = job[1]['vial_number'], \
                            sample_volume = self.config['redissolution']['solvent_volume']*ml, dilution = 1, solvent = solvent, water_ratio = 0, job = job[1]) 

                    num_redisol = 0 #reset counter
                    waiting_list = [] #reset waiting list
                    os.remove(processing_folder + job[0]) #remove job file
                    time.sleep(3)

            time.sleep(5)


    def measure_blank(self, fname, measurements, sample_position, sample_volume = 0.2*ml, solvent = 'ACN', **kwargs):


        self.selector = Valco_valve('visa://10.22.1.20/ASRL17::INSTR', dev_id = 0, mode = 3, position = 1)
        self.valve = Valco_valve('visa://10.22.1.20/ASRL4::INSTR', dev_id = 0, mode = 1, position = 'A')

        self.selector.move_to(1)

        self.log('washing selector to pump line with solvent {timestamp}')
        self.draw_and_dispense('ACN', 'valve', 0.3*ml, 500, 500)        

        # self.selector.move_to(6) #air

        self.log('Discard selector to pump volumes')
        # self.draw_and_dispense('valve', 'waste', 0.3*ml, 500, 500)
        self.draw_and_dispense('air', 'valve', 0.3*ml, 500, 500)  

        self.update_vial_status(sample_position, 'used')

        self.selector.move_to(sample_position)
        self.draw_and_dispense('ACN', 'valve', sample_volume + selector_to_pump, draw_velocity=500, dispense_velocity=500)
        time.sleep(3)

        self.do_measurements(fname, measurements = measurements, sample_position = sample_position, sample_volume = sample_volume, dilution = 1, solvent = solvent, water_ratio = 0, **kwargs) 


