from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time

class OptSigma_GSC02(VisaInstrument): 
    """  
    A Class control to the opti sigma's motorized stages by using GSC02.

    Parameters
    --------------------
    visa: str
        e.g. 'ASRL6::INSTR'
    res1 : float
        resolution of the stage (axis1)
    res2 : float
        resolution of the stage (axis2)
    verbose : bool
    """
    
    def __init__(self, visa, *args, res1 = 0.0025, res2 = 0.0025, verbose =False, **kwargs):

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
        
        self.verbose = verbose
        self.res1 = res1
        self.res2 = res2
        self.default_speed = {
            'range' : 2,                # 1: low speed # 2: high speed
            'axis1_min' :  500,         #pps (low:1-200, high:50-20000)
            'axis1_max' :  5000,        #pps (low:1-200, high:50-20000)
            'axis1_accel' : 200,        #ms 0-1000
            'axis2_min' :  500,         #pps (low:1-200, high:50-20000)
            'axis2_max' :  5000,        #pps (low:1-200, high:50-20000)
            'axis2_accel' : 200         #ms 0-1000
        } 

    def _PosToPulse(self, axis, pos):
        self._isValidAxis(axis, both = False)
        if axis == 1:
            return round(float(pos)/self.res1)
        elif axis == 2:
            return round(float(pos)/self.res2)     

    def _PulseToPos(self, axis, pulse):
        self._isValidAxis(axis, both = False)
        if axis == 1:
            return self.res1 * int(pulse) 
        elif axis == 2:
            return self.res2 * int(pulse) 
            

    def _isValidAxis(self, axis, both = True):
        if both == True:
            if axis not in [1,2,'W']:
                raise ValueError('axis should be 1, 2, or W')
        else:
            if axis not in [1,2]:
                raise ValueError('axis should be 1, 2')


    def _wait_ready(self):
        while True:
            if self.ask('!:') == 'R':
                break
            time.sleep(0.1)

    def set_speed(self, **kwargs):
        """  
        set moving speed of the stage

        Parameters
        --------------------
        params : 
             'range' : int (1: low speed, 2: high speed (2))
             'axis1_min' : int (low:1-200, high:50-20000 pps (500))
             'axis1_max' :  int (low:1-200, high:50-20000 pps (5000))
             'axis1_accel' : int (0-1000 ms (200))
             'axis2_min' :  int (low:1-200, high:50-20000 pps (500))
             'axis2_max' :  int (low:1-200, high:50-20000 pps (5000))
             'axis2_accel' : int (0-1000 ms (200))
        """
        setting = self.default_speed
        for key in kwargs.keys():
            setting[key] = kwargs[key]
        cmd = 'D:%sS%sF%sR%sS%sF%sR%s' %(setting['range'],
                                       setting['axis1_min'],
                                       setting['axis1_max'],
                                       setting['axis1_accel'],
                                       setting['axis2_min'],
                                       setting['axis2_max'],
                                       setting['axis2_accel']
                                       )
        self.write(cmd)


    def move_by(self, axis, step):
        """  
        Move stage to relative position  

        Parameters
        --------------------
        axis: int or str
            1, 2
        step: float
            moving length in the unit defined in the resolution
        """
        self._isValidAxis(axis, both = False)
        if step < 0:
            direction = '-'
        else : direction = '+'
        cmd = 'M:%s%sP%s' %(axis, direction, abs(self._PosToPulse(axis, step)))
        self.write(cmd)
        self.write('G')
        self._wait_ready()
        if self.verbose:
            self.get_position()


    def move_manual(self, axis, step):
        while True:
            command = input('input direction "c"(CW)/"r"(CCW) or position or new step "s+step" (current step %s)'%step)
            if command == 'c':
                self.move_by(axis, -step)
            elif command == 'r':
                self.move_by(axis, step)
            elif command[0] == 's':
                new_step = command.lstrip('s')
                if float(new_step) < 360:
                    step = float(new_step)
            elif command[0].isdigit() == True and float(command) < 360:
                self.move_to(axis, float(command))           
            elif command == 'q':
                break
            

    def move_by_both(self, step1, step2):
        """  
        Move both of the stage to relative position  

        Parameters
        --------------------
        step1: float
            moving length for axis1 in the unit defined in the resolution
        step2: float
            moving length for axis2 in the unit defined in the resolution
        """
        if step1 < 0:
            dir1 = '-'
        else : dir1 = '+'
        if step2 < 0:
            dir2 = '-'
        else : dir2 = '+'
        cmd = 'M:W%sP%s%sP%s' %(dir1, abs(self._PosToPulse(1, step1)), dir2, abs(self._PosToPulse(2, step2)))
        self.write(cmd)
        self.write('G')
        self._wait_ready()
        if self.verbose:
            self.get_position()


    def jog(self, axis):
        """  
        move stage until user stop(input s)

        Parameters
        --------------------
        axis: int or str
            1, 2
        """
        self._isValidAxis(axis, both = False) 
        direction = ''
        while True:
            direction = input('input "c"(CW), "r"(CCW) or "s"(Stop)')
            self.write('L:%s' %axis)
            self._wait_ready()
            if self.verbose:
                self.get_position()
            if direction == 's':
                break
            if direction == 'c':
                self.write('J:%s+' %axis)
                self.write('G')
            elif direction == 'r':
                self.write('J:%s-' %axis)
                self.write('G')
            
            

    # def manual_move(self, axis, default_step):
    #     """  
    #     Move stage manually with the key input

    #     Parameters
    #     --------------------
    #     axis: int or str
    #         1, 2
    #     default_step: float
    #         moving length in the unit defined in the resolution
    #     """
    #     self._isValidAxis(axis, both = False) 
    #     move = ''
    #     while True:
    #         move = input ('input "c"(CW), "r"(CCW), "s"(Stop) or pulse')
    #         step = self._PosToPulse(axis, default_step)

    #         if move == 's':
    #             break
    #         elif  move == 'c':
    #             direction = '+' 
    #         elif  move == 'r':
    #             direction = '-' 
    #         elif move < 0:
    #             step = abs(self._PosToPulse(axis, move))
    #             direction = '-'
    #         elif move >= 0:
    #             step = abs(self._PosToPulse(axis, move))
    #             direction = '+'
    #         cmd = 'M:%s%sP%s' %(axis, direction, abs(step))
    #         self.write(cmd)
    #         self.write('G')
    #         self._wait_ready()
    #         if self.verbose:
    #             self.get_position()


    def move_to(self, axis, position, cw = False):
        """  
        Move stage to absolute position  

        Parameters
        --------------------
        axis: int or str
            1, 2
        position: float
            destination in the unit defined in the resolution
        dir : bool
            stage always arrived to the new position from one direction if True
        """
        self._isValidAxis(axis, both = False) 

        if axis == 1:
            current_pos = self.get_position()[0]
        elif axis == 2:
            current_pos = self.get_position()[1]

        step = position - current_pos

        self.move_by(axis = axis, step = step)

        if cw and step < 0:
            self.move_by(axis=axis, step = - 1)
            self.move_by(axis=axis, step = 1)
    
    def move_to_both(self, pos1, pos2, cw = False):
        """  
        Move both of the stages to absolute position  

        Parameters
        --------------------
        pos1: float
            destination of the axis1 in the unit defined in the resolution
        pos2: float
            destination of the axis2 in the unit defined in the resolution
        """

        step1 = pos1 - self.get_position()[0]
        step2 = pos2 - self.get_position()[1]

        self.move_by_both(step1, step2) 

        if cw:
            if step1 < 0 and step2 < 0:
                self.move_by_both(-1, -1) 
                self.move_by_both(1, 1) 
            elif step1 < 0:
                self.move_by(axis = 1, step= -1) 
                self.move_by(axis = 1, step= 1)
            elif step2 < 0:
                self.move_by(axis = 2, step= -1) 
                self.move_by(axis = 2, step= 1)



    def move_home(self, axis, direction = '-'):
        """  
        Move stages to mechanical origin

        Parameters
        --------------------
        axis: int or str
            1, 2 or W(both)
        """
        self._isValidAxis(axis, both = True) 

        if axis == 1 or axis == 2:
            self.write('H:%s%s' %(axis, direction))
        elif axis == 'W':
            self.write('H:%s%s%s' %(axis, direction, direction))
        else:
            raise ValueError('axis should be 1, 2 or W')
        
        self._wait_ready()
        if self.verbose:
            self.get_position()


    def stop_move(self, axis):
        """  
        stop the stage movement

        Parameters
        --------------------
        axis: int or str
            1, 2 or W(both)
        """
        self._isValidAxis(axis, both = True) 
        self.write('L:%s' %axis)
        self._wait_ready()
        self.get_position()


    def get_position(self):
        """  
        get current position for each axis

        Returns
        --------------------
        position : list
            position of axis 1 and 2
        """
        pos_in_pulse = self.ask('Q').split(',')[0:2]
        position = []
        for i, pos in enumerate(pos_in_pulse):
            pos = pos.replace(' ', '')
            position.append(self._PulseToPos(i + 1, pos))
            # position.append(self._PulseToPos(i + 1, pos) % 360)
        # print('current position :\n axis1:%s, axis2 : %s'%(position[0], position[1]))
        return position


    def get_status(self):
        """  
        get current status for each axis

        Returns
        --------------------
        status : list
        """
        status = self.ask('Q').split(',')[2:]
        if status[0] == 'X':
            print('command or parameter errors')
        elif status[0] =='K':
            print('Command received normally')
        if status[1] == 'L':
            print('First axis stopped at limit switch') 
        elif status[1] == 'M':
            print('Second axis stopped at limit switch') 
        elif status[1] == 'W':
            print('Both axis stopped at limit switch') 
        elif status[1] == 'K':
            print('Normal Stop') 
        if status[2] == 'B':
            print('Busy')
        if status[2] == 'R':
            print('Ready')
        return status


    def set_zero(self, axis):
        """  
        set electrical origin for each axis

        Parameters
        --------------------
        axis: int or str
            1, 2 or W(both)
        """
        self._isValidAxis(axis, both = True)
        self.write('R:%s' %axis)
        self.get_position()


    def hold_motor(self, axis, hold = True):
        """  
        hold or free motor

        Parameters
        --------------------
        axis: int or str
            1, 2 or W(both)
        hold : bool
            hold motor if True
        """
        self._isValidAxis(axis, both = True) 
        self.write('C:%s%s' %(axis, hold))



