from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time


@add_set_get
class Harverd_11Elite(VisaInstrument):
    """  
    A class to control the Harverd syringe pump 11 Elite.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    """
    def __init__(self, visa, syringe_diameter, syringe_volume, *args, force = 50, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 115200,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r',
            'read_termination': '\r',
            'timeout': 1000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        self.diameter = syringe_diameter
        self.volume = syringe_volume

        print(self.ask('ver'))
        print(self.ask('address'))
        # print(self.ask('poll'))

        self.set_diameter(self.diameter)
        self.set_volume(self.volume)
        self.set_force(force)


    def ask(self, command):
        msg = self.manager.query(command)
        return msg.replace('\n','').replace(':', '')

    def set_force(self, force = 50):
        """
        set infusion force.
        see reference manual for recommended force for each syringe materials.
        (e.g. 50% for < 5ml, 100% for > 5ml for glass/plastic pumps)

        Parameters
        --------------------
        force : int
            infusion force 1-100
        """
        self.write('force %s' %force)
        time.sleep(0.1)
        #print(self.ask('force'))


    def set_diameter(self, diameter):
        """
        set diameter of syringe.

        Parameters
        --------------------
        diameter : float
            diameter of the syringe in mm
        """ 
        self.write('diameter %s mm' %diameter)
        time.sleep(0.1)
        #print(self.ask('diameter'))


    def set_volume(self, volume):
        """
        set volume of syringe.

        Parameters
        --------------------
        volume : float
            volume of the syringe in mL
        """ 
        self.write('svolume %s ml' %volume)
        time.sleep(0.1)
        #print(self.ask('svolume'))


    def set_infusion_rate(self, rate):
        """
        set infusion rate.

        Parameters
        --------------------
        rate : float
            rate of the syringe in uL/min
            (rate unit can be m, u, n, p/h, m, s)
        """          
        self.write('irate %s u/m' %rate)
        time.sleep(0.1)
        # self.ask('irate')
        # print('infusion rate : %s' %self.read())


    def set_withdrawal_rate(self, rate):
        """
        set withdrawal rate.

        Parameters
        --------------------
        rate : float
            rate of the syringe in uL/min
            (rate unit can be m, u, n, p/h, m, s)
        """          
        self.write('wrate %s u/m' %rate)
        # print(self.ask('wrate'))
        time.sleep(0.1)


    def set_target_volume(self, volume):
        """
        set target volume.

        Parameters
        --------------------
        volume : float
            target volume in uL.
        """          
        self.write('tvolume %s ul' %volume)


    def get_volume(self):
        """
        get infused/withdrawn volume.
        """        
        i_volume = self.ask('ivolume')
        w_volume = self.ask('wvolume')
        return [i_volume, w_volume]


    def clear_volume(self, direction):
        """
        clear the infused/withdrawn volume.

        Parameters
        --------------------
        direction : str
            "i", "w" or "both"
        """      
        if direction == 'i':
            self.write('civolume')
        elif direction == 'w':
            self.write('cwvolume')
        elif direction == 'both':
            self.write('cvolume')


    def clear_target_volume(self):
        """
        clear target volume.
        """         
        self.write('ctvolume')


    def set_target_time(self, time):
        """
        set target time.

        Parameters
        --------------------
        time : float
            target time in seconds.
        """          
        self.write('ttime %s' %time)


    def get_time(self):
        """
        get infused/withdrawn time.
        """        
        i_time = self.ask('itime')   
        w_time = self.ask('wtime')
        return [i_time, w_time]


    def clear_time(self, direction):
        """
        clear the infused/withdrawn time.

        Parameters
        --------------------
        direction : str
            "i", "w" or "both"
        """      
        if direction == 'i':
            self.write('citime')
        elif direction == 'w':
            self.write('cwtime')
        elif direction == 'both':
            self.write('ctime')


    def clear_target_time(self):
        """
        clear target time.
        """         
        self.write('cttime')


    def trigger_out_on(self):
        """
        set trigger out high
        """          
        self.write('output 1 high')


    def trigger_out_off(self):
        """
        set trigger out low
        """          
        self.write('output 1 low')


    def run(self, direction, rate, volume = None, duration = None, wait = False):
        """
        run the pump.

        Parameters
        --------------------
        direction : str
            "i" for infusion and "w" for withdrawal
        rate : float
             rate of the syringe in uL/min
        volume : float or None
            target volume in uL
        duration : float or None
            target time in seconds
        wait : Bool
            wait until pump stops if true
        """
        self.clear_time('both')
        self.clear_target_time()
        self.clear_volume('both')
        self.clear_target_volume()

        if volume:
            self.set_target_volume(volume)
        elif duration:
            self.set_target_time(duration)
        else:
            raise Exception('Harverd Elite 11:Please specify volume or duration')

        self.get_current_rate()

        if direction == 'i':
            self.set_infusion_rate(rate)
            self.write('irun')
            print('Harverd Elite 11: started to infuse %s uL at %s ul/min' %(volume, rate))
        elif direction == 'w':
            self.set_withdrawal_rate(rate)
            self.write('wrun')
            print('Harverd Elite 11: started to withdrawal %s uL at %s ul/min' %(volume, rate))
        else:
            raise Exception('Harverd 11 pump Error. Run direction should be "i" or "w".')

        if wait:
            flg = True
            while flg:
                if 'T*' in self.get_current_rate():
                    flg = False
                    print('Harverd Elite 11: Run finished')
                time.sleep(1)
                

    def stop(self):
        """
        Stop the pump.
        """
        self.write('stop')
        print('syringe stoped')


    def get_status(self):
        return self.ask('status')

    def get_current_rate(self):
        return self.ask('crate')

   