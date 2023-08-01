from pylab.instruments.pypowerusb.powerUSB import powerUSB as pUSB
import tkinter as tk
from tkinter import messagebox, ttk

class PowerManager():

    def __init__(self):
        
        self.PUSB = pUSB()
        self.num_device = self.PUSB.num_device

        self.dev_list = {
            'PM100D' : {'PowerUSB' : 0, 'port' : 1, 'measurement' : 'AbsPL'}, #power meter
            'QEPRO' : {'PowerUSB' : 0, 'port' : 2, 'measurement' : 'AbsPL'}, #spectrometer
            'SC20XR' : {'PowerUSB' : 0, 'port' : 3, 'measurement' : 'Evap'}, #shaker
            'DC4100' : {'PowerUSB' : 1, 'port' : 1, 'measurement' : 'AbsPL'}, #LED driver
            'KSC101' : {'PowerUSB' : 1, 'port' : 2, 'measurement' : 'AbsPL'}, #shutter controller
            'DH-MINI' : {'PowerUSB' : 1, 'port' : 3, 'measurement' : 'AbsPL'}, #D2/Halogen lamp
            'APD' : {'PowerUSB' : 2, 'port' : 1, 'measurement' : 'TE'}, #avalanche photo diode
            'PDL800D' : {'PowerUSB' : 2, 'port' : 2, 'measurement' : 'TE'}, #pico-sec laser driver
            'GSC02' : {'PowerUSB' : 2, 'port' : 3, 'measurement' : 'TA'} #sigma stage driver
        }

        self.dev_state = {
            '0' : 'off',
            '1' : 'on'
        }


    def list_device(self):

        dev_names = [key for key in self.dev_list.keys()]
        print(dev_names)
        return dev_names


    def select_device(self, device):

        if device not in self.dev_list:
            raise Exception ('PowerManager : %s is not in the device list' %device)
        else :
            dev = self.dev_list[device]

        self.PUSB.selectdevice(dev['PowerUSB'])

        return dev


    def read_state(self, device):
        """  
        Read power state of a specific device.

        Parameters
        --------------------
        device : str
            name of the device (should be in the device list)

        Returns
        --------------------
        state : int
            0(off) or 1(on)
        """
        dev = self.select_device(device)

        state = self.PUSB.readportstate()[dev['port'] - 1]

        print('PowerManager:: %s : %s ' %(device, self.dev_state[str(state)]))

        return state

    
    def set_state(self, device, state):
        """  
        Set power state of a specific device.

        Parameters
        --------------------
        device : str
            name of the device (should be in the device list)
        state : int (0 or 1) or str ('off' or 'on')
            0 : off, 1 : on
        """
        dev = self.select_device(device)

        if type(state) == str:
            for key, val in self.dev_state.items():
                if val == state.lower():
                    state = int(key)
                    break
        
        if state not in [0,1]:
            raise Exception('PowerManager : device state should be 0, 1 or "on", "off"')

        self.PUSB.setport(dev['port'], state)

        self.read_state(device)

        
    
