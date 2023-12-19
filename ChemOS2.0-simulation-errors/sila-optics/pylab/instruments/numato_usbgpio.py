from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const


@add_set_get
class Numato_Usbgpio(VisaInstrument):
    """  
    A class to control the Numato USB-GPIO.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    """
    def __init__(self, visa, *args, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 19200,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r',
            'read_termination': '\r',
            'timeout': 1000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        self.initialize()

    def initialize(self):
        self.ask('ver')
        print(self.read())


    def gpio_on(self, channel):
        """
        Turn on specific channel of gpio (5V).

        Parameters
        --------------------
        channel : int
            gpio channel to turn on  
        """
        self.ask('gpio set ' + str(channel))
   

    def gpio_off(self, channel):
        """
        Turn off specific channel of gpio (0V).

        Parameters
        --------------------
        channel : int
            gpio channel to turn off  
        """
        self.ask('gpio clear ' + str(channel))


    def gpio_read(self,channel):
        """
        read the state of specific channel of gpio.

        Parameters
        --------------------
        channel : int
            gpio channel to read

        Reterns
        --------------------
        state : int
            0(off) or 1(on)
        """
        self.ask('gpio read ' + str(channel))
        return self.read()

    def adc_read(self,channel):
        """
        read the adc input at specific channel of gpio.

        Parameters
        --------------------
        channel : int
            gpio channel to read 

        Reterns
        --------------------
        voltage : float
        """
        self.ask('adc read ' + str(channel))
        return 5*int(self.read())/1024




@add_set_get
class Numato_Usbrelay(VisaInstrument):
    """  
    A class to control the Numato USB-Relay.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    """
    def __init__(self, visa, *args, pos_count=None, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 19200,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\n\r',
            'read_termination': '\n\r',
            'timeout': 1000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        self.initialize()

    def initialize(self):
        self.ask('ver')
        print(self.read())

    def __del__(self):
        try:
            for i in range(4):
                self.relay_off(i)
            for i in range(4):
                self.relay_read(i)
            print('All relay off')
        except:
            pass
        

    def relay_on(self, channel):
        """
        Turn on specific channel of relay.

        Parameters
        --------------------
        channel : int
            relay channel to turn on  
        """
        self.ask('relay on ' + str(channel))
   

    def relay_off(self, channel):
        """
        Turn off specific channel of relay.

        Parameters
        --------------------
        channel : int
            relay channel to turn off  
        """
        self.ask('relay off ' + str(channel))


    def relay_read(self,channel):
        """
        read the state of specific channel of the relay.

        Parameters
        --------------------
        channel : int
            relay channel to read

        Reterns
        --------------------
        state : str
            off or on
        """
        self.ask('relay read ' + str(channel))
        state =  self.read()
        print('relay%s_' % channel +  state)
        return state


    def gpio_on(self, channel):
        """
        Turn on specific gpio channel (5V).

        Parameters
        --------------------
        channel : int
            gpio channel to turn on  
        """
        self.ask('gpio set ' + str(channel))
   

    def gpio_off(self, channel):
        """
        Turn off specific gpio channel(0V).

        Parameters
        --------------------
        channel : int
            gpio channel to turn off  
        """
        self.ask('gpio clear ' + str(channel))


    def gpio_read(self,channel):
        """
        read the state of specific channel of gpio.

        Parameters
        --------------------
        channel : int
            gpio channel to read

        Reterns
        --------------------
        state : int
            0(off) or 1(on)
        """
        self.ask('gpio read ' + str(channel))
        return self.read()
   
