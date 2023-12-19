from msl.loadlib import Server32
from ctypes import c_int, c_char, byref, create_string_buffer
import os

file_path = os.path.dirname(os.path.abspath(__file__))


class Myserver(Server32):
    def __init__(self, host, port, quiet, **kwargs):
        super(Myserver, self).__init__('%s/lib/PwrUSBDll.dll' %file_path, 'cdll', host, port, quiet)

    def InitPowerUSB(self):
        self.model = c_int()
        self.d = c_int()
        num_device = self.lib.InitPowerUSB(byref(self.model))
        return num_device, self.model.value

    def GetDeviceNumber(self):
        return self.lib.CheckStatusPowerUSB()

    def SelectDevice(self, dev_num):
        self.lib.SetCurrentPowerUSB(dev_num)

    def SetPortPowerUSB(self,port1, port2, port3):
        self.lib.SetPortPowerUSB(port1, port2, port3)

    def SetDefaultStatePowerUSB(self,port1,port2,port3):
        self.lib.SetDefaultStatePowerUSB(port1,port2,port3)

    def ClosePowerUSB(self):
        self.lib.ClosePowerUSB

    def ReadPortPowerUSB(self):
        port1, port2, port3 = c_int(), c_int(), c_int()
        self.lib.ReadPortStatePowerUSB(byref(port1), byref(port2), byref(port3))
        return port1.value, port2.value, port3.value

    def ReadDefaultPortState(self):
        port1, port2, port3 = c_int(), c_int(), c_int()
        self.lib.ReadDefaultPortStatePowerUSB(byref(port1), byref(port2), byref(port3))
        return port1.value, port2.value, port3.value

    def ReadCurrent(self):
        current = c_int()
        self.lib.ReadCurrentPowerUSB(byref(current))
        return current.value

