import ctypes as ct
import os,sys
from ctypes.util import find_library
from ctypes import c_long, c_int, byref
import time
from .thorlabs_APT import ThorlabsAPT


def list_available_device():
    """  
    get serial numbers of the available devices

    Returns
    --------------------
    serial_number : list of int
    """
    APT = ThorlabsAPT(serial_number = None)
    serial_number = APT.list_available_devices()
    print(serial_number)


# class ThorlabsK10CR1(ThorlabsAPT):
class ThorlabsK10CR1():
    """  
    This is the python module to control to the thorlabs rotation stage 
    K10CR1 by using ThorlabsAPT.

    Parameters
    --------------------
    serial_number: int
        serial number of the stage
    sethome : bool
        move to home when initializing the device if True
    """
    def __init__(self, serial_number, sethome = True):
        # super(ThorlabsK10CR1,self).__init__(serial_number)
        self.apt = ThorlabsAPT(serial_number)

        self.apt.enable()
        self.apt.set_move_home_parameters(2, 4, 8, 0)
        self.apt.set_velocity_parameters(0,8,10)
        if sethome:
            self.move_home()
            self.move_to(0)
            print('done')
       

    def _wait_move(self):
        while self.apt.is_in_motion():
            time.sleep(0.2)


    def move_home(self):
        """  
        move to home position
        """
        self.apt.move_home()
        self._wait_move()


    def move_by(self, value):
        """  
        move stage by using the relative position

        Parameters
        --------------------
        value: float
            angle in degrees
        """
        self.apt.move_by(value)
        self._wait_move()


    def move_to(self, value):
        """  
        move stage by using absolute position

        Parameters
        --------------------
        value: float
            angle in degrees
        """
        self.apt.move_to(value)
        self._wait_move()


    def set_velosity(self, max_vel, accel):
        """  
        set velosity of the stage

        Parameters
        --------------------
        max_val: float
            angle/seconds
        accel: float
            angle/seconds2
        """
        self.apt.set_velocity_parameters(0, max_vel, accel)


    def get_position(self):
        """  
        get current position

        Returns
        --------------------
        position : float
        """
        return self.apt.position
