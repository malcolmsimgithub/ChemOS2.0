from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time

class Valco_valve(VisaInstrument): 
    """  
    A Class control to the valco's injector and flow selector valves.
    When first time use valve, long press home button of remote actuator, select "2. Interface Setup", select "2. USB", select "2. 9600"

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
    
    def __init__(self, visa, *args, dev_id, mode, position = 'A', varbose = False,  **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r',
            'read_termination': '\r',
            'timeout': 5000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        self.varbose = varbose
        self.dev_id = dev_id
        self.write('ID%s' %self.dev_id)
        self.mode = mode
        if self.mode == 3:
            self.offset = int(self.ask_offset())
        self._set_response_mode(0)
        self._set_mode(mode)
        self.move_to(position)


    def _ask(self, command):
        self.write(command)
        return self.manager.read()[len(command):]


    def _set_mode(self, mode):
        self.manager.query('%sAM%s' %(self.dev_id, mode))


    def _set_response_mode(self, mode):
        """  
        set the response mode of the valve

        Parameters
        --------------------
        mode: int
            0: no response 
            1: basic response string
            2: extended response string
        """
        if self.mode == 1:
            self.write('%sIFM%s' %(self.dev_id, mode))
        if self.mode == 3:
            self.manager.query('%sIFM%s' %(self.dev_id, mode))
        #print(self.read())
#    def get_status(self):
#        print(self._ask('%sSTAT' %self.dev_id))
#        print(self.read(0))


    def move_to(self, position, waiting = 1):
        """  
        Move valve to the specified position

        Parameters
        --------------------
        position: str, or int
            target position('A" or 'B' for mode = 1, 2 and 1 to number of position for mode = 3)
        waiting : float
            waiting time after the command
        """
        #print(position)
        if position == 'A':
            self.write('%sCW' %self.dev_id)
        elif position == 'B':
            self.write('%sCC' %self.dev_id)
        else:
            self.write('%sGO%s' %(self.dev_id, position + self.offset -1))
        time.sleep(waiting)
        if self.varbose:
            print('valve %s moved to %s' %(self.dev_id, self.ask_position()))


    def move_home(self, waiting = 1):
        """  
        Move valve to home position

        Parameters
        --------------------
        waiting : float
            waiting time after the command
        """
        self.write('%sHM' %self.dev_id)
        time.sleep(waiting)

        if self.varbose:
            print('valve %s moved to %s' %(self.dev_id, self.ask_position()))

    def ask_position(self):
        """  
        ask current valve position

        Returns
        --------------------
        position : str
        """
        return self._ask('%sCP' %self.dev_id)


    def ask_offset(self):
        """  
        ask current offset value for valve position

        Returns
        --------------------
        position : str
        """
        return self._ask('%sSO' %self.dev_id)

    def get_count(self):
        """  
        get the moving count of the valve

        Returns
        --------------------
        count : int
        """
        return self._ask('%sCNT' %self.dev_id)

    def find_stop(self):
        self.write('%sLRN' %self.dev_id)

    def get_moving_time(self):
        """  
        get the time taken for the last movement of the valve

        Returns
        --------------------
        time(ms) : float
        """
        return self._ask('%sTM' %self.dev_id)

    def close_device(self):
        self.manager.close()
        print('valve "%s" closed' %self.dev_id)




