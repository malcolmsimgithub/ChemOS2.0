from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time

class Coherent_OBIS(VisaInstrument): 
    """  
    A Class control to the valco's injector and flow selector valves.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    dev_id : int (0~9)
    mode : int
        type of the valve (1:2position valves, 3:multiposition valves)
    position : str
        initial position
    """
    
    def __init__(self, visa, *args, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r\n',
            'read_termination': '\r\n',
            'timeout': 5000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)
        
        self.max_power = float(self.ask('SOURce:POWer:LIMit:HIGH?'))
        self.wavelength = self.get_wavelength()
        print('OBIS %s nm: baseplate temparature %s' %(self.wavelength, self.get_baseplate_temperature()))
        #self.laser_off()


    def _ask(self, command):
        self.ask(command)
        return self.manager.read()


    def laser_on(self, power = None, power_rel = None):
        if power == None:
            power = self.max_power * power_rel

        self.ask('SOURce:POWer:LEVel:IMMediate:AMPLitude %s' %power)
        self.ask('SOURce:AM:STATe ON')

        print('OBIS %s nm : on (%s mW, BaseTemp: %s)' %(self.wavelength, float(self.get_power())*1000, self.get_baseplate_temperature()))

    def laser_off(self):
        self.ask('SOURce:AM:STATe OFF')
        print('OBIS %s nm : off' %(self.wavelength))

    def set_power(self, power = None, power_rel = None):
        if power == None:
            power = self.max_power * power_rel
        self.ask('SOURce:POWer:LEVel:IMMediate:AMPLitude %s' %power)

    def get_power(self):
        return self._ask('SOURce:POWer:LEVel:IMMediate:AMPLitude?')

    def get_status(self):
        return self._ask('SOURce:AM:STATe?')


    def get_wavelength(self):
        return self._ask('SYSTem:INFormation:WAVelength?')

    def get_baseplate_temperature(self):
        return self._ask('SOURce:TEMPerature:BASeplate?')