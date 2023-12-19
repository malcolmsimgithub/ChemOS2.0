"""
Module: instruments

    This module consists classes for managing common instruments in Soft Semiconductor Group. Please refer to the class for detailed description. The basic use is declaring the instrument with visa, and calling its method to use, similar to using the instrument directly.

    Instruments currently supported:
        Magnet
            - Magnet Controller, Lakeshore 642
            - Methods:
        Lockin
            - Lock In Amplifier, SR830
        MagnetMeter
            - Magnet Meter, Lakeshore 425
        Monochromator
        Keithley
            - Power source, multimeter, Keithley 2400
        FunctionGenerator
            - Function Generator, HP33120A
        Timer
            - System timing instrument
        HP
        CryoController
            - CroCen 32B
        Oscilloscope
            - Tektronix TDS 3054C



Manager : Tony Wu
Email   : tonyw@mit.edu
License : Open source, under MIT License
Author  : Tony Wu
Credits : Tony Wu, Nick Thompson, Dan Congreve
Updates :
    - 2016/03/18 updated list_instruments(), new pyvisa raises new errors or COM ports
    - 2015/10/16 updated HP, fixed timeout issues on measurements. OLED EQE measurement program should work out of box.
    - 2015/9/16 updated to new pyvisa, better handling with serial instruments

TODO    : Test Keithley.set_measure_range()
"""


import visa
import pyvisa.constants as pvc
import time, math





### Functions
def list_instruments():
    "Lists the available visa resources and its corresponding instruments."
    visa_lists = visa.ResourceManager().list_resources()
    for v in visa_lists:
        instrument_model = None
        try:
            instrument_v = Instrument(v)
            instrument_model = instrument_v.test()[:-2]
            print "%-18s%-5s%s" % (v, ">>>", instrument_model)
        except (visa.VisaIOError):
            print "%-18s%-5s%s" % (v, ">>>", "Cannot get model information.")



def get_instruments():
    "Get the available visa resources and its corresponding ID."
    visa_lists = visa.ResourceManager().list_resources()
    instrument_list = []
    for v in visa_lists:
        instrument_model = None
        try:
            instrument_v = Instrument(v)
            instrument_model = instrument_v.test()[:-2]
            instrument_list.append({
                'visa': v,
                'model': instrument_model.split(',')
            })
        except (visa.VisaIOError):
            instrument_list.append({
                'visa': v,
                'model': [None]
            })
    return instrument_list
















### Base Instrument Manager Class
class Instrument:
    """
    Base Instrument Manager Class for different instruments
    """

    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        self.__manager = None
        if visa_string:
            self.set_visa(visa_string, **kwargs)
        else:
            self.__manager = None

    # Load Instrument by Visa
    def set_visa(self, visa_string, **kwargs):
        "Sets the instrument's visa"
        if visa_string:
            self.__manager = visa.ResourceManager().open_resource(visa_string, **kwargs)

    def test(self):
        "Test the model of __manager"
        return str(self.__manager.ask('*IDN?'))

    # IO
    def write(self, input_string):
        "Sends commands to the instrument"
        if input_string:
            self.__manager.write(input_string)

    def read(self):
        "Reads instrument values"
        return self.__manager.read()

    def ask(self, query_string):
        "write(query_string) and returns read()"
        if query_string:
            return self.__manager.ask(query_string)
        return

    # Methods that must be defined in inherited classes
    def read_values(self):
        "Default read_values() in base instrument class"
        return

















### Specific Instrument Controllers with high level control

class Magnet(Instrument):
    """
    Magnet Controller, Lakeshore 642
        Manual:
            http://www.lakeshore.com/Documents/642_Manual.pdf

        Example:
            >>> my_magnet = instruments.Magnet('GPIB1::12::INSTR')
            >>> my_magnet.set_current(30)
            This will run the magnet on 30 amps (~0.35 Tesla) for 12 loops, approximately 5 minutes.
    """
    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        self.mode_dict = {
            'OFF': 0,
            'ON': 1,
            'AUTO': 2,
            'DISABLED': 3
        }

        Instrument.__init__(self, visa_string, **kwargs)


    def __del__(self):
        "Sets currnt to 0 A, protects Lakeshore 642 if the program exits improperly"
        self.set_current(0)
        self.set_internal_water('AUTO')

    # Defining the inherited method read_values()
    def read_values(self):
        "Read the current on the magnet"
        return float(self.ask('RDGI?'))

    def set_current(self, current):
        "Set magnet current (amps)"
        if current is not None:
            self.write('SETI %.4f' % current)


    def get_internal_water(self):
        "Get internal water mode"
        mode = int(self.ask('INTWTR?'))
        modes = ['OFF', 'ON', 'AUTO', 'DISABLED']
        return modes[mode]

    def set_internal_water(self, mode):
        "Set internal water mode, OFF, ON, AUTO, DISABLED"
        if isinstance(mode, str):
            mode = self.mode_dict[mode.upper()]
        self.write('INTWTR %d' % mode)

    def cooling(self, seconds = 10):
        mode = self.get_internal_water()
        time.sleep(1)
        self.set_internal_water('on')
        for t in range(seconds, 0, -1):
            print 'Cooling Lakeshore ... {:2d}s\r'.format(t),
            time.sleep(1)
        print ''
        self.set_internal_water(mode)





class Magnet6674A(Instrument):
    """
    Magnet Controller, Agilent 6674A
        Manual:
            http://literature.cdn.keysight.com/litweb/pdf/5964-8269.pdf?id=1000002298-1:epsg:man

        Example:
            >>> my_magnet = instruments.Magnet6674A('GPIB0::4::INSTR')
            >>> my_magnet.set_current(30)
            This will run the magnet on 30 amps (~0.35 Tesla) for 12 loops, approximately 5 minutes.
    """
    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        Instrument.__init__(self, visa_string, **kwargs)
        self.set_voltage(15)
        self.set_current(0)
        self.on()


    def __del__(self):
        "Sets current to 0 A, protects PSU if the program exits improperly"
        self.off()
        self.set_voltage(0)
        self.set_current(0)
        

    # Defining the inherited method read_values()
    def read_values(self):
        "Read the current on the magnet"
        return float(self.ask('MEASURE:CURRENT?'))

    def set_current(self, current):
        "Set magnet current (amps)"
        if current is not None:
            self.write('CURRENT %.2f' % current)

    def set_voltage(self, voltage):
        "Set magnet voltage (volts)"
        if voltage is not None:
            self.write('VOLTAGE %.2f' % voltage)

    def on(self):
        self.write('OUTPUT ON')
    def off(self):
        self.write('OUTPUT OFF')












class Lockin(Instrument):
    """
    Lock In Amplifier SR830
        Manual:
            http://www.thinksrs.com/downloads/PDFs/Manuals/SR830m.pdf

        Example:
            >>> my_lockin = instruments.Lockin('GPIB1::3::INSTR')
            >>> my_lockin.set_sensitivity('Auto')
            >>> my_lockin.read_values()
    """

    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        # Lockin Variables
        self.__wait_time = 1.0 # 1 second
        self.__sens_high_threshold = 0.9
        self.__sens_low_threshold  = 0.36
        self.__read_list = ['X', 'Y', 'R', 'T']

        Instrument.__init__(self, visa_string, **kwargs)

    # Defining the inherited method read_values()
    def read_values(self, read_list = None):
        """
        Read lockin values (in list). read_list accepts list of integers. Options include: [1, 2, 3, 4] representing [X, Y, R, Theta]. The function also supports string options ['X', 'Y', 'R', 'T']. Additional options include [5, 6, 7, 8, 9, 10, 11] = Aux In 1~4 and Reference Frequency, CH1, CH2. Reads default values if read_list is not presented.
        """
        options_dict = {'x': 1, 'y': 2, 'r': 3, 't': 4}
        if read_list is None:
            read_list = self.__read_list

        if isinstance(read_list, list):
            # Scan throught read_list and change strings to corresponding value
            for i in range(0, len(read_list)):
                read_list[i] = self.__check_dict(match_dict = options_dict, value = read_list[i], default = 3)

            # Determine which command to use with different read length
            if len(read_list) == 1:
                return [float(self.ask('OUTP? %d;' % read_list[0]))]
            elif len(read_list) <= 6:
                lockin_values = self.ask( 'SNAP? ' + ','.join(str(values) for values in read_list) +';' )
                return map(float, lockin_values.split(','))
            else:
                print "Length for read_list is too long, please input lists less than 6 elements."
                return

        else:
            print "read_list must be of data type list"
            return

    def set_read_list(self, read_list):
        "Set the default values to read by read_values()"
        if read_list is not None:
            self.__read_list = read_list

    # Setting Methods
    def set_coupling_mode(self, mode):
        "Set lockin coupling mode, AC (default) or DC"
        if mode:
            if mode.lower() == 'dc':
                self.write('ICPL 1;')
            else: # other ac
                self.write('ICPL 0;')

    def set_query_mode(self, mode):
        "Set lockin measuring mode: A, A-B, I (low), I (high), corresponding to mode = 0, 1, 2, 3. Allows string options: 'A', 'A-B', 'Ilow', 'Ihigh'"
        mode_dict = {'a': 0, 'ab': 1, 'ilow': 2, 'ihigh': 3}
        mode = self.__check_dict(match_dict = mode_dict, value = mode)
        if mode is not None:
            self.write('ISRC %d;' % mode)

    def set_query_noise(self, noise):
        "Set query noise: noise = 0, 1, 2 (High Reserve, Normal, Low Noise). Allows string options: 'High', 'Normal', 'Low'"
        noise_dict = {'high': 0, 'normal': 1, 'low': 2}
        noise = self.__check_dict(match_dict = noise_dict, value = noise)
        if noise is not None:
            self.write('RMOD %d;' % noise)

    def set_display(self, channel, display = 0):
        "Set the front panel display. Channel 1, 2; display 0 (channels) or 1 (X/Y)"
        # DDEF?
        if channel is not None:
            self.write('FPOP %d,%d;' % (channel, display))

    def set_sensitivity(self, sensitivity):
        """
        Set sensitivity values from 0~26 corresponding to 2fA ~ 1uA or 2nV ~ 1V (every 2, 5, 10 incrementals). Allows string input including unit, no space (eg. '2uV', '1pA', '5nA'). Use 'Auto' for best sensitivity by searching automatically. The automatic method is inspired by Nick Thompson's code.
        """
        curr_dict = self.__make_dict(['2fa', '5fa', '10fa', '20fa', '50fa', '100fa', '200fa', '500fa', '1pa', '2pa', '5pa', '10pa', '20pa', '50pa', '100pa', '200pa', '500pa', '1na', '2na', '5na', '10na', '20na', '50na', '100na', '200na', '500na', '1ua'])
        volt_dict = self.__make_dict(['2nv', '5nv', '10nv', '20nv', '50nv', '100nv', '200nv', '500nv', '1uv', '2uv', '5uv', '10uv', '20uv', '50uv', '100uv', '200uv', '500uv', '1mv', '2mv', '5mv', '10mv', '20mv', '50mv', '100mv', '200mv', '500mv', '1v'])
        curr_volt_dict = dict(curr_dict.items() + volt_dict.items() + [('auto', 'auto')])
        sensitivity = self.__check_dict(match_dict = curr_volt_dict, value = sensitivity)
        if sensitivity is not None:
            # auto mode
            while sensitivity == 'auto':
                r_value = self.read_values(['R'])[0]
                base_value = 1e-9 if (int(self.ask('ISRC?;')) <= 1) else 1e-15
                cur_sens_index = int(self.ask('SENS?;'))
                cur_sens_value = self.__index_to_value( cur_sens_index, [2, 5, 10], base_value)
                if r_value > (cur_sens_value * self.__sens_high_threshold) and cur_sens_index < 26:
                    self.write('SENS %d;' % (cur_sens_index+1) )
                elif r_value < (cur_sens_value * self.__sens_low_threshold) and cur_sens_index > 0:
                    self.write('SENS %d;' % (cur_sens_index-1) )
                elif r_value > (cur_sens_value * self.__sens_high_threshold) and cur_sens_index >= 26:
                    print( 'Signal Saturated, please lower intensity!' )
                    sensitivity = 26
                else:
                    sensitivity = cur_sens_index
            self.write('SENS %d;' % sensitivity)

    def set_automatic_sensitivity_thresholds(self, high, low):
        "Setting the thresholds (high, low) for automatic search sensitivity."
        self.__sens_high_threshold = high
        self.__sens_low_threshold  = low

    # Time related settings
    def set_time_constant(self, time_constant):
        """
        Set time constant with options 0~13 (corresponding to 10us~30ks for every 1, 3). Also allows string input including unit, no space (eg. '100us', '10ms', '3s'). Default to 100ms if string not found. Automatically sets to the best wait time for read value.
        """
        tc_dict = self.__make_dict(['10us', '30us', '100us', '300us', '1ms', '3ms', '10ms', '30ms', '100ms', '300ms', '1s', '3s', '10s', '30s', '100s', '300s', '1ks', '3ks', '10ks', '30ks'])
        time_constant = self.__check_dict(match_dict = tc_dict, value = time_constant)
        if time_constant is not None:
            self.write('OFLT %d;' % time_constant)
            self.set_wait_time()

    def set_low_pass_filter(self, low_pass_filter):
        """
        Set low pass filter 6, 12, 18, 24 dB/oct (0, 1, 2, 3). Accepts both 0~3 or strings with unit 'dB', no space (eg. Lockin.set_filter('24dB')). Automatically sets to the best wait time for read value.Automatically set the best wait time.
        """
        db_dict = {'6db': 0, '12db': 1, '18db': 2, '24db': 3}
        low_pass_filter = self.__check_dict(match_dict = db_dict, value = low_pass_filter)
        if low_pass_filter is not None:
            self.write('OFSL %d;' % low_pass_filter)
            self.set_wait_time()

    def set_wait_time(self, wait_time = None):
        """
        Sets the wait time for reading values based upon filter and time constants of the lockin. Also allows manual input, note that changing time constant and filter will automatically change the wait time.
        """
        if wait_time is not None:
            self.__wait_time = wait_time
        else:
            wait_time_multiple = {'0': 5.0, '1': 7.0, '2': 9.0, '3': 10.0}
            low_pass_filter    = int(self.ask('OFSL?;'))
            time_constant      = self.__index_to_value( int(self.ask('OFLT?;')), [1, 3], 1e-5)
            self.__wait_time   = wait_time_multiple[str(low_pass_filter)] * time_constant


    def set_query_rate(self, rate):
        "Set the data sampling rate from 62.5mHz (0) to 512Hz (13) or Trigger (14). Also allows string input including unit, no space (eg. '125mHz', '256Hz', 'trigger')."
        rate_dict = self.__make_dict(['62.5mhz', '125mhz', '250mhz', '500mhz', '1hz', '2hz', '4hz', '8hz', '16hz', '32hz', '64hz', '128hz', '256hz', '512hz', 'trigger'])
        rate = self.__check_dict(match_dict = rate_dict, value = rate)
        if rate is not None:
            self.write('SRAT %d;' % rate)

    # Buffer related Settings
    def set_buffer(self, mode):
        "Set buffer mode: Shot (0) or Loop (1). Loop mode will erase the oldest buffer to read new values. Shot will stop reading if buffer is full."
        if mode is not None:
            self.write('SEND %d;' % mode)

    def clear_buffer(self):
        "Clears Lockin buffers."
        self.write('REST;')

    def reset_buffer(self):
        "Reset the lockin data buffers and set buffer to loop mode."
        self.set_buffer(1)
        self.clear_buffer()

    # Reading Lockin
    def get_wait_time(self):
        "Returns the optimized time to wait before reading Lock In. It is strongly suggested to always wait at least this amount of time between consecutive lockin measurement."
        return self.__wait_time

    def is_overload(self):
        "Query if output, returns True if overloading."
        self.write('*CLS')
        return bool(int(self.ask('LIAS?2;')))


    ## Private Utility Functions for Lockin Methods
    def __make_dict(self, string_list, start_int = 0):
        "Returns dictionary to that match string list to integers."
        return_dict = dict()
        list_length = len(string_list)
        for i in range(0, list_length):
            return_dict[ string_list[i] ] = i + start_int
        return return_dict

    def __check_dict(self, match_dict, value, default = None):
        "Check if the value is string and convert to value by match_dict."
        if isinstance(value, str):
            return match_dict.get(value.lower(), default)
        elif isinstance(value, int):
            return value
        else:
            return default


    def __index_to_value(self, index, value_structures, base_value, base_index = 0):
        "Calculates the index to values, used for lockin values"
        index = index - base_index
        struct_len = len(value_structures)
        # Gets the order and structure to use
        order_of_value = index / struct_len
        struct_index= index % struct_len

        return base_value * float(value_structures[struct_index]) * 10.0**order_of_value

        # # Read with single value request
        # if isinstance( read_list, int):
        #     return float(self.ask('OUTP? %d;' % read_list))
        # elif isinstance( read_list, str):
        #     read_list = self.__check_dict(match_dict = options_dict, value = read_list, default = 3)
        #     return float(self.ask('OUTP? %d;' % read_list))

        # # Read with multiple values request
        # elif isinstance(read_list, list) and len(read_list) > 1 and len(read_list) <= 6:
        #     # scan through if there are strings and convert them to corresponding intergers
        #     for i in range(0, len(read_list)):
        #         read_list[i] = self.__check_dict(match_dict = options_dict, value = read_list[i], default = 3)
        #     lockin_values = self.ask( 'SNAP? ' + ','.join(str(values) for values in read_list) +';' )
        #     return map(float, lockin_values.split(','))
        # else:
        #     return 0

















#  Inherits from Instrument
class MagnetMeter(Instrument):
    """
    Magnet Flux Meter, Lakeshore 425
        Manual: http://www.lakeshore.com/Documents/425Manual.pdf

        Example:
            >>> my_magnet_meter = MagnetMeter('COM7')
            >>> my_magnet_meter.read_values()

        Serial port configuration:
            baud_rate = 57600
            data_bits = 7
            stop_bits = 1
            parity = 1 (odd)
            term_chars = '\\r\\n'

        Troubleshooting:
            Please close all other programs accessing the magnet meter, as serial ports does not support multiple connections.
    """

    # Defining the inherited method read_values()
    # Needs to fix the error in the first read.
    def read_values(self):
        "Reads Magnetic Field"
        magnetic_field = self.ask('RDGFIELD?')
        # Clearing Buffer
        if magnetic_field == '':
            magnetic_field = self.ask('RDGFIELD?')
        return float(magnetic_field)
        # return self.ask('RDGFIELD?')

    # Overwrites constructor due to serial connection
    def __init__(self, visa_string = None):
        """
        Constructor with optional visa input
        """
        self.magnet_meter_settings = {
            'baud_rate': 57600,
            'data_bits': 7,
            'stop_bits': pvc.StopBits.one,
            'parity'   : pvc.Parity.odd,
            'end_input': pvc.SerialTermination.termination_char
        }
        Instrument.__init__(self, visa_string, **self.magnet_meter_settings)


    # Load Instrument by Visa
    def configure_connection(self):
        "(Obsolete) Sets the instrument's visa"
        # Parameters for only Lakeshorer 642
        self._Instrument__manager.timeout = 10
        self._Instrument__manager.baud_rate = 57600
        self._Instrument__manager.data_bits = 7
        self._Instrument__manager.stop_bits = 1
        self._Instrument__manager.parity = 1
        self._Instrument__manager.term_chars = '\r\n'

        # self._Instrument__manager.baud_rate = 9600
        pass





#  Inherits from Instrument
class Monochromator(Instrument):
    """
    Monochromator, Oriel Cornertone 130
        Manual: https://www.newport.com/medias/sys_master/images/images/h8f/hb0/8797226860574/Oriel-Cornerstone-130-User-Manual-RevA.pdf

        Example:
            >>> mc = Monochromator('COM2')
            >>> mc.read_values()

        Serial port configuration:
            baud_rate = 9600
            data_bits = 8
            stop_bits = 1
            parity = none
            term_chars = '\\r\\n'

        Troubleshooting:
            Please close all other programs accessing the monochromator.
    """



    # Overwrites constructor due to serial connection
    def __init__(self, visa_string = None):
        """
        Constructor with optional visa input
        """
        self.monochromator_settings = {
            'baud_rate': 9600,
            'data_bits': 8,
            'stop_bits': pvc.StopBits.one,
            'parity'   : pvc.Parity.none,
            'end_input': pvc.SerialTermination.termination_char,
        }
        Instrument.__init__(self, visa_string, **self.monochromator_settings)

    # Overrides write ask
    def write(self, input_string):
        "Sends commands to the instrument"
        if input_string:
            echo = Instrument.ask(self, input_string)

    def ask(self, query_string):
        "write(query_string) and returns read()"
        if query_string:
            echo = Instrument.ask(self, query_string)
            return self.read()
        return

    # methods
    def read_values(self):
        "Reads wavelength "
        return self.ask('WAVE?')

    def go_wave(self, wavelength):
        self.write('GOWAVE {:3.2f}'.format(wavelength))

    def units(self, units=None):
        "Query units with no options. Set units with options include NM (nanometers), UM (micrometers), WN (wavenumbers)"
        if units == None:
            return self.ask('UNITS?')
        else:
            self.write('UNITS {}'.format(units))

    # shutter
    def shutter(self, t=None):
        "Query shutter status if no arguments. Set shutter status with open or close."
        if t == None:
            status = self.ask('SHUTTER?')
            if status[0] == 'O':
                return 'open'
            else:
                return 'close'
        elif t[0].upper() == 'O':
            self.write('SHUTTER O')
        elif t[0].upper() == 'C':
            self.write('SHUTTER C')
        else:
            print 'Error shutter status: use O or C.'
        return


















class Keithley(Instrument):
    """
    Keithley Model 2400
        Manual:

        Example:
            >>> my_keithley = Keithley('GPIB1::24::INSTR')
            >>> my_keithley.set_measure_range_auto('On')
            >>> my_keithley.set_source_value(0.1)
            >>> my_keithley.on()
            >>> my_keithley.read_values()
    """
    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        self.__mode = 'VOLT'
        Instrument.__init__(self, visa_string, **kwargs)

    # Defining the inherited method read_values()
    def read_values(self, read_list = None):
        "Reads values from Keithley, default is VOLT,CURR,RES,TIME,STAT, change format with set_measure_values()."
        if read_list is not None:
            self.set_read_list(read_list)
        return map(float, self.ask('READ?').split(','))

    # Turn Keithley on off
    def on(self):
        "Turn Keithley on"
        self.write('OUTP:STAT ON')
    def off(self):
        "Turn Keithley off"
        self.write('OUTP:STAT OFF')

    # Source Settings
    def set_source_mode(self, mode):
        """
        Set source mode, CURRent or VOLTage. This also initialize the keithley to measure fixed values and syncing measurement and display range. You can use write() to overwrite some predefined settings manually.
        """
        if mode is None:
            return

        abv_mode = mode[0:4].upper()
        if abv_mode == 'CURR' or abv_mode == 'VOLT':
            self.__mode = abv_mode
        else:
            self.__mode = 'VOLT'

        # Setting the mode (and setting the measure mode to opposite)
        self.write('SOUR:FUNC %s' % self.__mode)
        measure_mode = 'CURR' if self.__mode == 'VOLT' else 'VOLT'
        self.write('SENS:FUNC "%s"' % measure_mode)
        self.write('SOUR:%s:MODE FIX' % self.__mode)
        self.write('SENS:%s:PROT:RSYN ON' % self.__mode)


    def set_source_range_auto(self, onoff):
        "Turning on/off for auto source range"
        if onoff is not None:
            self.write('SOUR:%s:RANG:AUTO %s' % (self.__mode, onoff.upper()))

    def set_source_range(self, source_range):
        "Set source range manually"
        if source_range is not None:
            if source_range.upper() == 'AUTO':
                self.set_source_range_auto('ON')
            else:
                self.write('SOUR:%s:RANG %f' % (self.__mode, source_range))

    def set_source_value(self, value):
        "Set source output value"
        if value is not None:
            self.write('SOUR:%(mode)s:IMM %(value)f' % {'mode': self.__mode, 'value': value})

    # Measurement Settings
    def set_measure_range_auto(self, onoff):
        "Turning on/off for auto measuring range"
        measure_mode = 'CURR' if self.__mode == 'VOLT' else 'VOLT'
        if onoff is not None:
            self.write('SENS:%s:RANG:AUTO %s' % (measure_mode, onoff.upper()))

    def set_measure_range(self, measure_range):
        "Set measurement range"
        if measure_range is not None:
            if measure_range.upper() == 'AUTO':
                self.set_measure_range_auto('ON')
            else:
                self.write('SENS:%s:RANG %f' % (measure_mode, measure_range))

    def set_read_list(self, read_list):
        """"
        Set the values to read. read_list takes list of strings. Default is all values, read_list = ['VOLT', 'CURR', 'RES', 'TIME', 'STAT']. Calling set_measure_values() without read_list reset to default.
        """
        if read_list is None:
            return

        # Check for valid input
        for r in read_list:
            if r.upper() not in ['VOLT', 'CURR', 'RES', 'TIME', 'STAT']:
                print "Invalid value to read: %s" % r
                return
        self.write('FORM:ELEM %s' % ','.join([r.upper() for r in read_list]) )
























class FunctionGenerator(Instrument):
    """
    Function Generator: HP33120A
        Manual:
            https://www.ini.uzh.ch/~ppyk/BasicsOfInstrumentation/HP_33120A_Function_Generator.pdf

        Example:
            >>> fg = FunctionGenerator('GPIB1::10::INSTR')
            >>> fg.apply_pulse_signal(frequency = 213, voltage = 5)
    """

    # Defining the inherited method read_values()
    def read_values(self):
        "Reads current settings from function generator."
        shape       = str(self.ask('FUNC:SHAP?'))[:-1]
        frequency   = float(self.ask('FREQ?'))
        vpp         = float(self.ask('VOLT?'))
        offset      = float(self.ask('VOLT:OFFS?'))
        return [shape, frequency, vpp, offset]

    # High Level Commands
    def apply_signal(self, shape, frequency, vpp, offset):
        """
        Apply Signal with function shape (shape), frequency (freq), peak to peak voltage (vpp), signal offset (offset).
        Options for shapes includs: sinusoid (SIN), square (SQU), triangle (TRI), RAMP, noise (NOIS), DC.
        Frequency unit is Hz. Voltage units are V.
        """
        if shape and frequency and vpp and offset:
            self.write("APPL:%(shape)s %(frequency)f, %(vpp)f, %(offset)f" % {'shape': shape.upper(), 'frequency': frequency, 'vpp': vpp, 'offset': offset})

    def apply_pulse_signal(self, frequency, voltage):
        "Apply pulse signal with frequency and peak voltage."
        if frequency and voltage:
            self.apply_signal(shape = 'SIN', frequency = frequency, vpp = voltage, offset = float(voltage)/2)

    def turn_off(self):
        "Turn function generator off (to minimum values)."
        self.set_offset(0.0)
        self.set_vpp(0.05)

    # Low Level Commands
    def set_shape(self, shape):
        "Set signal shape, options: sinusoid (SIN), square (SQU), triangle (TRI), RAMP, noise (NOIS), DC."
        if shape:
            self.write('FUNC:SHAP %s' % shape.upper())
    def set_frequency(self, frequency):
        "Set frequency, unit Hz"
        if frequency:
            self.write('FREQ %f' % frequency)
    def set_vpp(self, vpp):
        "Set peak to peak voltage, units V."
        if vpp:
            self.write('VOLT %f' % vpp)
    def set_offset(self, offset):
        "Set signal offset, units V."
        if offset:
            self.write('VOLT:OFFS %f' % offset)

    def clear_status(self):
        "Clear Error Status"
        self.write('*CLS')
















class Timer():
    """
    Timer
        Example:
        >>> t = Timer()
        >>> t.read_values()
        >>> t.reset_time()
    """

    def __init__(self):
        self.__start_time = time.time()

    def read_values(self):
        "Return time elapsed."
        return time.time() - self.__start_time

    def reset_time(self):
        "Resets the timer time."
        self.__start_time = time.time()





##### ADDED in Bawendi's Lab #####
class Keithley2602(Instrument):
    """
    Keithley Model 2602
        Manual:
            http://www.ece.uprm.edu/~etclab/resources/equipment/keithley2612/2600series_referencemanual.pdf

        Example:
            >>> my_keithley = Keithley2602('GPIB0::26::INSTR')
            >>> my_keithley.set_measure_range_auto('ON')
            >>> my_keithley.set_source_value(0.1)
            >>> my_keithley.on()
    """
    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        self.__mode = 'VOLTS'
        Instrument.__init__(self, visa_string, **kwargs)

    # Defining the inherited method read_values()
    def read_values(self, read_list = None):
        pass

    # Turn Keithley on off
    def on(self):
        "Turn Keithley on"
        self.write('smua.source.output = smua.OUTPUT_ON')
    def off(self):
        "Turn Keithley off"
        self.write('smua.source.output = smua.OUTPUT_OFF')

    # Source Settings
    def set_source_mode(self, mode):
        """
        Set source mode, VOLTS and AMPS.
        """

        if mode.upper() in ['VOLTS', 'AMPS']:
            self.__mode = mode.upper()
        else:
            raise ValueError( 'Wrong mode. Please use VOLTS or AMPS')

        # Setting the mode (and setting the measure mode to opposite)
        self.write('smua.source.func = smua.OUTPUT_DC%s' % self.__mode)


    def set_source_range_auto(self, onoff):
        "Turning on/off for auto source range"
        if onoff.upper() in ['ON', 'OFF']:
            self.write('smua.source.autorangev = smua.AUTORANGE_' + onoff.upper())
            self.write('smua.source.autorangei = smua.AUTORANGE_' + onoff.upper())
        else:
            raise ValueError ('Wrong onoff value')

    def set_source_range(self, source_range):
        "Set source range"
        if source_range == 'AUTO':
            self.set_source_range_auto('ON')
        else:
            if self.__mode == 'VOLTS':
                self.write('smua.source.rangev = {}'.format(value))
            else:
                self.write('smua.source.rangei = {}'.format(value))

    def set_source_value(self, value):
        "Set source output value"
        if self.__mode == 'VOLTS':
            self.write('smua.source.levelv = {}'.format(value))
        else:
            self.write('smua.source.leveli = {}'.format(value))

    # Measurement Settings
    def set_measure_range_auto(self, onoff):
        "Turning on/off for auto measuring range"
        if onoff.upper() in ['ON', 'OFF']:
            self.write('smua.measure.autorangev = smua.AUTORANGE_' + onoff.upper())
            self.write('smua.measure.autorangei = smua.AUTORANGE_' + onoff.upper())
        else:
            raise ValueError ('Wrong onoff value')

    def set_measure_range(self, measure_range):
        "Set measurement range"
        if measure_range == 'AUTO':
            self.set_measure_range_auto('ON')
        else:
            if self.__mode == 'VOLTS':
                self.write('smua.measure.rangev = {}'.format(value))
            else:
                self.write('smua.measure.rangei = {}'.format(value))

    def set_read_list(self, read_list):
        """"
        Set the values to read. read_list takes list of strings. Default is all values, read_list = ['VOLT', 'CURR', 'RES', 'TIME', 'STAT']. Calling set_measure_values() without read_list reset to default.
        """
        pass



class SorensenDLM(Instrument):
    """
    SorensenDLM
        Manual:

        Example:
            >>> sorensen = Keithley2602('GPIB0::26::INSTR')
            >>> sorensen.set_measure_range_auto('ON')
            >>> sorensen.set_source_value(0.1)
            >>> sorensen.on()
    """
    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        Instrument.__init__(self, visa_string, **kwargs)

    # Defining the inherited method read_values()
    def read_values(self):
        return self.measure('current')

    # output on/off
    def on(self):
        self.write('OUTPUT:STATE ON')
    def off(self):
        self.write('OUTPUT:STATE OFF')

    # set compliance
    def set_limit(self, voltage = None, current = None):
        if voltage != None:
            self.write('SOURCE:VOLTAGE:LIMIT {}'.format(voltage))
        if current != None:
            self.write('SOURCE:CURRENT:LIMIT {}'.format(current))

    # set voltage and current
    def output(self, voltage = None, current = None):
        if voltage != None:
            self.write('SOURCE:VOLTAGE {}'.format(voltage))
        if current != None:
            self.write('SOURCE:CURRENT {}'.format(current))

    # measure
    def get(self, mode):
        "Measures Voltage or Current."
        if mode.upper()[0] == 'V':
            return float(self.ask('MEASURE:VOLTAGE?'))
        elif mode.upper()[0] == 'C':
            return float(self.ask('MEASURE:CURRENT?'))
        else:
            raise ValueError('mode should be VOLTAGE or CURRENT')

class BrushedMotor():
    def __init__(self, visa, library_folder = r'C:\Program Files\Thorlabs\Kinesis', mm_index = 25./857600.):
        # load libraries
        os.environ['PATH'] += os.pathsep + library_folder
        self._kdll = ct.CDLL('Thorlabs.MotionControl.KCube.DCServo.dll')
        self._mm_index = mm_index

        # check if serial is in device list
        self._kdll.TLI_BuildDeviceList()
        if self._kdll.TLI_GetDeviceListSize() < 1:
            raise ValueError('No Device found.')
        devices_string = ct.c_char_p(' '*100)
        self._kdll.TLI_GetDeviceListByTypeExt(devices_string, ct.c_int(100), ct.c_int(27))
        if str(visa) not in devices_string.value:
            raise ValueError('Serial number {} not found in device list'.format(visa))

        # create words
        self._mtype = ct.wintypes.WORD()
        self._mid   = ct.wintypes.WORD()
        self._mdata = ct.wintypes.DWORD()

        # open port
        self._serial = ct.c_char_p(str(visa))
        if self._kdll.CC_Open(self._serial) != 0:
            raise ValueError('Cannot open {}'.format(visa))
        self._kdll.CC_StartPolling(self._serial, ct.c_int(100)) # set polling to 100ms

        # Clear some weird message
        self._kdll.CC_ClearMessageQueue(self._serial)
        self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))
        self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))

    def to_home(self):
        self._kdll.CC_ClearMessageQueue(self._serial)
        self._kdll.CC_Home(self._serial)

        self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))
        while (self._mtype.value != 2) or (self._mid.value != 0):
            self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))

        return self

    def to_position(self, position):
        'Go to position in mm'
        self._kdll.CC_ClearMessageQueue(self._serial)
        self._kdll.CC_MoveToPosition(self._serial, ct.c_int(int(position/self._mm_index)))

        self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))
        while (self._mtype.value != 2) or (self._mid.value != 1):
            self._kdll.CC_WaitForMessage(self._serial, ct.byref(self._mtype), ct.byref(self._mid), ct.byref(self._mdata))

        return self

    def read_values(self):
        'Get the position in mm'
        return self.get_position()
    def get_position(self):
        'Get the position in mm'
        return self._kdll.CC_GetPosition(self._serial) * self._mm_index

    def __del__(self):
        "Delete object"
        self._kdll.CC_StartPolling(self._serial)
        self._kdll.CC_Close(self._serial)


##### END #####




class HP(Instrument):
    """
    HP: HP 4156C
        Manual:
            http://cp.literature.agilent.com/litweb/pdf/04156-90050.pdf
            http://wiki.epfl.ch/carplat/documents/4156_SCPI_Command_Reference.pdf

        Example:
        >>> hp = HP('GPIB0::19::INSTR')
        >>> hp.set_I(ch=1, i=-2E-4).read_values()
    """
    icomp = 1e-2 # 10mA
    vcomp = 25

    def __init__(self, visa_string = None):
        """
        Constructor with optional visa input
        """
        if visa_string:
            self.set_visa(visa_string, timeout=25000)
        else:
            self.__manager = None

    def reset(self):
        "Resets HP. For detailed format and modes, it is written in page 1-140 ~ 1-105."
        self.write('*RST')
        self.write('US')
        self.write('FMT 2, 1')
        return self

    def set_integration_time(self, time):
        "Sets measure integration time to short (S), medium (M), long (L)."
        ti = 'SML'.find(time[0].upper()) + 1
        ti = 1 if ti is -1 else ti
        self.write('SLI ' + str(ti))
        return self

    def set_channels(self, channels):
        "Set the HP to measure the channels, supports single channel or list of channels"
        # check if channels are list
        if not isinstance(channels, list):
            channels = [channels]
        # set modes
        self.set_mode(channels)
        # close unused channels
        off_channels = [c for c in [1,2,3,4] if not c in channels]
        self.turn_off_channel(off_channels)
        # open channels
        self.turn_on_channel(channels)
        return self


    def set_V(self, ch, v, icomp = None):
        "Setting bias voltage for measurement. Current compliance uses default value (%e) if not entered." % self.icomp
        if icomp == None:
            icomp = self.icomp
        self.write('DV %d,0,%s,%s' % (ch, str(v), str(icomp)) )
        return self

    def set_I(self, ch, i, vcomp = None):
        "Setting bias current for measurement. Voltage compliance uses default value (%f) if not entered." % self.vcomp
        if vcomp == None:
            vcomp = self.vcomp
        self.write('DI %d,0,%s,%s' % (ch, str(i), str(vcomp)) )
        return self


    def read_values(self):
        "Read and return values from HP"
        self.write('XE')
        string_values = self.ask('RMD?')
        values = map(float, string_values.split(','))
        return values

    def save_data(self, f):
        variables = self.ask(':PAGE:DISP:LIST?').replace('\n', '').split(',')
        # write to file
        out_str = ' '.join(['{}'] * (len(variables)+1)) + '\n'
        units = []
        data = []
        for var in variables:
            units += [self.ask(":DATA:UNIT? '{}'".format(var)).replace('\n', '')]
            data += [self.ask(":DATA? '{}'".format(var)).replace('\n', '').split(',')]
        
        # write header
        headers = ['No.'] + ['{}({})'.format(v, u) for v, u in zip(variables, units)]
        f.write(out_str.format(*headers))

        # write data
        for i in range(0, len(data[0])):
            f.write( out_str.format(*([i+1]+[d[i] for d in data])) )
        return



    # Functions for more detailed options, use set_channel for easier and packaged usages

    def set_mode(self, channels):
        "Setting channels to be measured"
        if not isinstance(channels, list):
            channels = [channels]
        self.write('MM 1,' + ','.join(map(str, channels)))
        return self

    def turn_on_channel(self, channels):
        "Turn the specified channels on. May use list of integers or single channel number. Turns off other channels."
        if not isinstance(channels, list):
            channels = [channels]
        self.write('CN ' + ','.join(map(str, channels)))
        return self

    def turn_off_channel(self, channels):
        "Turn off the specified channels on. May use list of integers or single channel number. Turns off other channels."
        if not isinstance(channels, list):
            channels = [channels]
        self.write('CL ' + ','.join(map(str, channels)))
        return self

    def set_Icomp(self, icomp):
        "Change default current compliance"
        self.icomp = icomp
        return self

    def set_Vcomp(self, vcomp):
        "Change default voltage compliance"
        self.vcomp = vcomp
        return self








class CryoController(Instrument):
    """
    CryoCon 32B
        Manual:
            http://www.cryocon.com/Model32/M32UM.pdf

        Example:
            >>> cc = CryoController('GPIB1::11::INSTR')
            >>> cc.
    """

    def __init__(self, visa_string = None, **kwargs):
        """
        Constructor with optional visa input
        """
        self.channel = 'A'
        self.loop = 1
        Instrument.__init__(self, visa_string, **kwargs)

    def read_values(self):
        "Reads Temperature"
        return float(self.ask('INPUT {}:TEMP?'.format(self.channel)))

    def set_channel(self, channel = None):
        "Sets channel. A or B."
        if channel:
            self.channel = channel
        return self

    def set_loop(self, loop = None):
        "Sets loop. 1 or 2."
        if loop:
            self.loop = loop
        return self

    def set_temperature(self, temperature = None):
        "Set Loop Temperature"
        if temperature >= 0:
            self.write('LOOP {}:SETPT {}'.format(self.loop, temperature))
        return self

    def get_set_temperature(self):
        "Get Set Loop Temperature"
        return float(self.ask('LOOP {}:SETPT?'.format(self.loop))[:-1])

    def set_pid(self, p=None, i=None, d=None):
        "Set PID Loop parameter. set_pid(p=5, i=5, d=0)"
        if p >= 0:
            self.write('LOOP {}:PGAIN {}'.format(self.loop, p))
        if i >= 0:
            self.write('LOOP {}:IGAIN {}'.format(self.loop, i))
        if d >= 0:
            self.write('LOOP {}:DGAIN {}'.format(self.loop, d))
        return self

    def set_range(self, lrange = None):
        "Set Loop Range"
        if lrange.upper() not in ['HI', 'MID', 'LOW']:
            raise ValueError('Please type in HI, MID or LOW')
        if lrange.upper():
            self.write('LOOP {}:RANGE {}'.format(self.loop, lrange.upper()))
        return self

    def get_range(self):
        "Get Loop Range"
        return self.ask('LOOP {}:RANGE?'.format(self.loop))






class Oscilloscope(Instrument):
    """
    Tektronix TDS 3054C
        Manual:
            http://ics-web.sns.ornl.gov/dht/TEK-3000-series-programming-manual.pdf

        Example:
            >>> oc = Oscilloscope('GPIB1::1::INSTR')
    """

    def __init__ (self, *arg, **kwargs):
        "Overide initialization"
        Instrument.__init__(self, *arg, **kwargs)
        self.set_read_settings()

    def set_mode(self, mode):
        "Set Acquiring mode: SAMPLE, AVERAGE, PEAK, ENVELOPE"
        if mode:
            self.write('ACQUIRE:MODE {}'.format(mode))
        return self

    def set_averages(self, num):
        "Set Average numbers"
        if num:
            self.write('ACQUIRE:NUMAVG {:d}'.format(num))
        return self

    def set_read_settings(self, source = 1, start = 1, stop = 100000, width=2, encoding = 'ASCII'):
        "Set Read Settings."
        if source:
            self.write('DATA:SOURCE CH{:d}'.format(source))
        if start:
            self.write('DATA:START {:d}'.format(start))
        if stop:
            self.write('DATA:STOP {:d}'.format(stop))
        if width:
            self.write('DATA:WIDTH {:d}'.format(width))
        if encoding:
            self.write('DATA:ENC {:s}'.format(encoding))
        return self

    def get_channels(self):
        "Get channels that are on"
        select_str = self.ask('SEL?').split(';')[0:4]
        channels = []
        for i in range(0, 4):
            if int(select_str[i]) == 1:
                channels += [i+1]
        return channels

    def set_y(self, source = 1, volts = None, pos = None):
        "Set Y scope"
        if volts != None:
            self.write('CH{}:VOLTS {}'.format(source, volts))
        if pos != None:
            self.write('CH{}:POS {}'.format(source, pos))
        return self

    def get_y(self, source = 1):
        "Get Y scope"
        volts = self.ask('CH{}:VOLTS?'.format(source))
        pos = self.ask('CH{}:POS?'.format(source))
        return volts, pos

    def start(self):
        "Start the scope continuously"
        self.write('ACQUIRE:STOPAFTER RUNSTOP')
        self.write('ACQUIRE:STATE RUN')
        return self

    def stop(self):
        "Stop the scope"
        self.write('ACQUIRE:STATE STOP')
        return self


    def run(self, block = True):
        "Take a single acquisition based on mode and average numbers. By default, will block program if not finished."
        self.write('ACQUIRE:STOPAFTER SEQ')
        self.write('ACQUIRE:STATE RUN')
        while block and (self.ask('ACQUIRE:STATE?')[0] == '1'):
            time.sleep(0.1)

    def read_curve(self):
        "Reads data"
        data = self.ask('CURVE?').split(',')
        ymult = float(self.ask('WFMP:YMULT?'))
        yzero = float(self.ask('WFMP:YZERO?'))
        data = [float(c)*ymult+yzero for c in data]
        xincr = float(self.ask('WFMP:XINCR?'))
        xzero = float(self.ask('WFMP:XZERO?'))
        t = [xzero+i*xincr for i in range(0, len(data))]
        return [t, data]





class PulseGenerator(Instrument):
    """
    HP 8114A
        Manual:
            http://cp.literature.agilent.com/litweb/pdf/5980-1213E.pdf

        Example:
            >>> pg = PulseGenerator('GPIB1::13::INSTR')
    """

    def on(self):
        "Turn on"
        self.write('OUTPUT:STATE ON')
        return self

    def off(self):
        "Turn off"
        self.write('OUTPUT:STATE OFF')
        return self

    def set_frequency(self, frequency):
        "Set frequency, in Hz"
        self.write('SOURCE:FREQUENCY {}'.format(frequency))

    def set_duty_cycle(self, duty_cycle):
        "Set duty cycle, available inputs from 1~100 (%)"
        self.write('SOURCE:PULSE:DCYCLE {}'.format(duty_cycle))

    def set_voltage(self, low, high):
        "Set low, high voltages"
        self.write('SOURCE:VOLTAGE:LOW {}'.format(low))
        self.write('SOURCE:VOLTAGE:HIGH {}'.format(high))










### Test

def main(argv=None):
    return # success

if __name__ == '__main__':
    status = main()
    sys.exit(status)
