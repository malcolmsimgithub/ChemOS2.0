from msl.loadlib import Client64
import time
import os.path

class powerUSB(Client64):
    """ 
    This is a Python module to control the powerUSB(http://www.pwrusb.com/).
    """

    def __init__(self):
        # super(powerUSB, self).__init__(module32 = 'C:/Users/MatterLab/Dropbox/PythonLab/pylab/instruments/pypowerusb/my_server')
        module_folder = os.path.dirname(__file__)
        super(powerUSB, self).__init__(module32 = os.path.join(module_folder, 'my_server'))

        self.num_device, model = self._initPUSB()
        print('PowerUSB : number of PowerUSB connected %s' %self.num_device)

    def _initPUSB(self):
        """  
        initialize the connection to the powerusb.
        """
        return self.request32('InitPowerUSB')

    def getdevicenumber(self):
        self.num_device = self.request32('GetDeviceNumber')
    
    def selectdevice(self, dev_num):
        self.request32('SelectDevice', dev_num)

    def readcurrent(self):
        current = self.request32('ReadCurrent')
        return current

    def readportstate(self):
        """  
        Read states of the ports

        Reterns
        --------------------
        states : list
            list of the port state (0: off, 1: on)
        """
        port1, port2, port3 = self.request32('ReadPortPowerUSB')
        return [port1, port2, port3]

    def readdefaultstate(self):
        """  
        Read default states of the ports

        Reterns
        --------------------
        states : list
            list of the port state (0: off, 1: on)
        """
        port1, port2, port3 = self.request32('ReadDefaultPortState')
        return [port1, port2, port3]    


    def setport(self, port, state = 0):
        """  
        Set states of a specific port.

        Parameters
        --------------------
        port : int
            port to be set (1-3)
        state : int
            0 : off, 1 : on
        """
        new_state = self.readportstate()
        for i in range(3):
            if i == port - 1:
                new_state[i] = state
        
        self.request32('SetPortPowerUSB', new_state[0], new_state[1], new_state[2])

        time.sleep(1)


    def setport_all(self, port1, port2, port3):
        """  
        Set states of all the ports.

        Parameters
        --------------------
        port1 : int
            0: off, 1: on
        port2 : int
            0: off, 1: on
        port3 : int
            0: off, 1: on
        """
        self.request32('SetPortPowerUSB', port1, port2, port3)

        time.sleep(1)


    def setdefault(self, port1, port2, port3):
        """  
        Set default states of the ports.

        Parameters
        --------------------
        port1 : int
            0: off, 1: on
        port2 : int
            0: off, 1: on
        port3 : int
            0: off, 1: on
        """
        self.request32('SetDefaultStatePowerUSB', port1, port2, port3)


    def close(self):
        """  
        close connection to the powerusb.
        """
        self.request32('ClosePowerUSB')
