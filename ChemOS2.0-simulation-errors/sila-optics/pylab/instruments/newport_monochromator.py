from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time

@add_set_get
class CS210(VisaInstrument):
    """  
    This class controls the Newport CS210 monochromator

    Parameters
    --------------------
    visa : str
        visa address of CS210
    """
    
    def __init__(self, visa, *args, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'end_input': pv_const.SerialTermination.termination_char,
#            'write_termination' : '\n',
#            'read_termination' : '\r\n',
            'timeout': 5000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

    # Overrides write ask
    def send(self, input_string):
        "Sends commands to the instrument"
        if input_string:
            echo = VisaInstrument.ask(self,input_string)
            print(echo)

    def query(self, query_string):
        "write(query_string) and returns read()"
        if query_string:
            echo = VisaInstrument.ask(self,query_string)
            return self.read()
        return


    # methods
    def read_values(self):
        """  
        Reads wavelength

        Returns
        --------------------
        wavelength : float
        """
        return self.query('wave?')


    def go_wave(self, wavelength):
        """  
        move grating to specific wavelength

        Parameters
        --------------------
        wavelength : float
        """
        self.send('GOWAVE {:3.2f}'.format(wavelength))


    def units(self, units=None):
        """  
        Query units with no options. 
        Set units with options include NM (nanometers), UM (micrometers), WN (wavenumbers)

        Parameters
        --------------------
        units : str, or None
            NM (nanometers), UM (micrometers), WN (wavenumbers)

        Returns
        --------------------
        units : str
        """
        if units == None:
            return self.query('UNITS?')
        else:
            self.send('UNITS {}'.format(units))


    #grating
    def grating(self, number=None):
        """  
        Query active grating number if no arguments.
        Set grating with options.

        Parameters
        --------------------
        number : int
            grating number to be set

        Returns
        --------------------
        number : int
            active grating number 
        """
        if number == None:
            return self.query('GRAT?')
        else:
            self.send('GRAT ' + str(number) )
            time.sleep(17)

    # shutter
    def shutter(self, t=None):
        """  
        Query shutter status if no arguments.
        Set shutter status with open or close.

        Parameters
        --------------------
        t : str
            'open' or 'close'

        Returns
        --------------------
        t : str
            'open' or 'close'
        """
        
        if t == None:
            status = self.query('SHUTTER?')
            if status[0] == 'O':
                return 'open'
            else:
                return 'close'
        elif t[0].upper() == 'O':
            self.send('SHUTTER O')
        elif t[0].upper() == 'C':
            self.send('SHUTTER C')
        else:
            print ('Error shutter status: use O or C.')
        return



##self.initialize(pos_count)
#
#
#    def device_info(self):
#        return self.ask('*idn?')
#
#    def initialize(self, pos_count):
#        self.write('sensors=0')
#        self.write('speed=1')
#        if pos_count:
#            self.set_position_count(pos_count)
#
#    # ---- override methods ----
#    def write(self, command):
#        if command:
#            self.manager.write(command)
#            return self.read()
#
#    def set_position_count(self, pos_count):
#        self.write('pcount={}'.format(pos_count))
#    def get_position_count(self):
#        return int(self.ask('pcount?'))
#
#    @rangemethod(1, 'get_position_count', int)
#    def set_position(self, pos):
#        self.write('pos={}'.format(pos))
#    def get_position(self):
#        return int(self.ask('pos?'))
    
