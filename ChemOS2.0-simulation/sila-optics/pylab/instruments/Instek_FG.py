from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const


@add_set_get
class Instek_FG(VisaInstrument):
    """  
    A class to control the Instek funtion generator.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    """
    def __init__(self, visa, *args, **kwargs):

        #rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\n',
            'read_termination': '\n',
            'timeout': 1000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)
        
        self.clb = 0.084

        # self._initialize()
        
    def _initialize(self):
        self.write('*idn?')
        print(self.read())


    def sin_out(self, frequency, Amplitude, offset):
        """  
        Output the sin wave with calibrated amplitude and offset. 

        Parameters
        --------------------
        frequency : float
            output frequency (Hz)
        Amplitude : float
            output amplitude (V)
        offset : float
            offset of the wave (V) 
        """
        self.write('SOURce1:APPL:SIN %s,%s,%s' %(frequency, (Amplitude/(1+self.clb)), (offset/(1+self.clb/2))))


    def square_out(self, frequency, Amplitude, offset, duty):
        """  
        Output the square wave with calibrated amplitude and offset. 

        Parameters
        --------------------
        frequency : float
            output frequency (Hz)
        Amplitude : float
            output amplitude (V)
        offset : float
            offset of the wave (V)
        duty : float
            duty ratio of the output wave (%) 
        """
        self.write('SOUR1:FUNC SQU')
        self.write('SOUR1:FREQ %s' % frequency)
        self.write('SOUR1:AMPL %s' % (Amplitude/(1+self.clb)))
        self.write('SOUR1:DCO %s' % (offset/(1+self.clb)))
        self.write('SOUR1:SQU:DCYC %s' % duty)
        self.write('OUTP ON')


    def square_out_arb(self, frequency, Amplitude, offset, duty, delay = 0):
        """  
        Output the square wave with calibrated amplitude and offset
        by using arbitrary waveform.
        Square wave with the duty less than 1% can be generated with this function.
        (mininum duty : 0.05%) 

        Parameters
        --------------------
        frequency : float
            output frequency (Hz)
        Amplitude : float
            output amplitude (V)
        offset : float
            offset of the wave (V)
        duty : float
            duty ratio of the output wave (%) 
        delay : float
            delay of the signal (0-1)
        """
        self.arbit_set(self._createStepFunction(4000, duty*0.01, delay = delay))
        self.arbit_out(frequency, Amplitude *2, offset - Amplitude/2)


    def ramp_out(self, frequency, Amplitude, offset):
        """  
        Output the ramp wave with amplitude and offset. 

        Parameters
        --------------------
        frequency : float
            output frequency (Hz)
        Amplitude : float
            output amplitude (V)
        offset : float
            offset of the wave (V)
        """
        self.write('SOURce1:APPL:RAMP %s,%s,%s' %(frequency, Amplitude, offset))


    def _createStepFunction(self, datalength, duty, delay):
        """  
        create step function for arbitrary output.

        Parameters
        --------------------
        datalength: int
            length of the data for the waveform
        duty : float
            duty ratio (0-1)
        delay : float 
            delay of the signal (0-1)
    
        Reterns
        --------------------
        wave_form : list
        """
        if duty < 0 or 1 < duty:
            raise ValueError('values of the duty should be {} to {}'.format(0, 1))
        elif delay < 0 or 1 < delay:
            raise ValueError('values of the delay should be {} to {}'.format(0, 1))

        waveform = [0 for i in range(datalength)]
        for i in range(int(datalength*duty)):
            waveform[i+int(datalength*delay)] = 1

        return waveform


    def arbit_set(self, waveform):
        """  
        set user defined arbitrary waveform to the volatile memory of the FG.

        Parameters
        --------------------
        waveform: list
            value should be -1 to 1. 
            data length should be 2 to 4096.
        """
        if len(waveform) < 2 or len(waveform) > 4096 :
            raise ValueError('waveform length should be {} to {}'.format(2, 4096))

        if min(waveform) < -1 or max(waveform) > 1:
            raise ValueError('values of the waveform should be {} to {}'.format(-1, 1))

        # cmd = 'DATA:DAC VOLATILE, 0, %s' %waveform
        cmd = ('DATA:DAC VOLATILE, 0, %s' %[round(d*511) for d in waveform]).replace('[', '').replace(']','')
        self.write(cmd)


    def arbit_out(self, frequency, Amplitude, offset):
        """  
        Output the user defined wave with calibrated amplitude and offset. 

        Parameters
        --------------------
        frequency : float
            output frequency (Hz)
        Amplitude : float
            output amplitude (V)
        offset : float
            offset of the wave (V)
        """
        self.write('SOUR1:FREQ %s' % frequency)
        self.write('SOUR1:AMPL %s' % (Amplitude/(1+self.clb)))
        self.write('SOUR1:DCO %s' % (offset/(1+self.clb)))
        self.write('SOUR1:FUNC USER')
        self.write('OUTP ON')


    def off(self):
        """  
        turn off the output signals
        """
        self.write('OUTP OFF')

    
    def save_setting(self, mem_num):
        """  
        save current device setting to the memory.

        Parameters
        --------------------
        mem_num : int
            memory number where the data is saved. (0-9)
        """
        self.write('*SAV %s' %mem_num)


    def save_arbitrary_save(self, mem_num):
        """  
        save user defined wave to the memory.

        Parameters
        --------------------
        mem_num : int
            memory number where the data is saved. (10-19)
        """
        self.write('*SAV %s' %mem_num)  


    def load_setting(self, mem_num):
        """  
        load device setting pr user defined wave.

        Parameters
        --------------------
        mem_num : int
            memory number where the data is saved. (0-19)
        """
        self.write('*RCL %s' %mem_num)  