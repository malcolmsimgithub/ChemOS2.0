#### using Kinesis
from .. import CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
from ctypes import cdll, create_string_buffer
import ctypes as ct
# from ctypes.wintypes import WORD, DWORD
import os, time

__all__ = ['ThorlabsK10CR1']

class DeviceInfo(ct.Structure):
    # _fields_ = [
    #     ('typeID', DWORD),
    #     ('description', ct.c_char*65),
    #     ('serialNo', ct.c_char*9),
    #     ('PID', DWORD),
    #     ('isKnownType', ct.c_bool),
    #     ('maxChannels', ct.c_short),
    # ]
    _fields_ = [
        ('typeID', ct.c_ulong),
        ('description', ct.c_char*65),
        ('serialNo', ct.c_char*9),
        ('PID', ct.c_ulong),
        ('isKnownType', ct.c_bool),
        ('maxChannels', ct.c_short),
    ]

module_path = os.path.dirname(os.path.abspath(__file__))
dll_path = os.path.join(module_path, 'dll')

@add_set_get
class ThorlabsK10CR1(object):
    def __init__(self, visa, set_home = True, dllpath=dll_path, 
            dllname=r'Thorlabs.MotionControl.IntegratedStepperMotors.dll'):

        os.environ['PATH'] = dllpath + os.pathsep + os.environ['PATH']
        self.apt = cdll.LoadLibrary(dllname)

        # check visa num
        if self.apt.TLI_BuildDeviceList() != 0:
            raise ValueError('Something wrong with device')
        data = create_string_buffer(100)
        self.apt.TLI_GetDeviceListByTypeExt(data, 100, 55)
        if visa in str(data.value):
            self.serial = ct.c_char_p(str.encode(visa))
            self.excmd('ISC_Open')
        else:
            self.serial = ct.c_char_p(str.encode(''))
            raise ValueError('Device Not Found')

        self.set_vel_params(10, 10)
        if set_home:
            self.home()


    def __del__(self):
        self.excmd('ISC_Close')

    def wait_message(self, msg=(2, 0)):
        msg_type = ct.c_ushort(10)
        msg_id = ct.c_ushort(10)
        msg_data = ct.c_ulong()
        while msg_type.value != msg[0] or msg_id.value != msg[1]:
            self.excmd('ISC_WaitForMessage', ct.byref(msg_type), ct.byref(msg_id), ct.byref(msg_data))
            time.sleep(0.1)
        return msg_data.value

    def excmd(self, command, *args):
        resp = getattr(self.apt, command)(self.serial, *args)
        time.sleep(0.1)
        return resp

    def device_info(self):
        info = DeviceInfo()
        self.excmd('TLI_GetDeviceInfo', ct.byref(info))
        return '{}, serial: {}'.format(info.description.decode(), info.serialNo.decode())

    def clear_message(self):
        self.excmd('ISC_ClearMessageQueue')

    def home(self):
        self.clear_message()
        self.excmd('ISC_Home')
        self.wait_message(msg=(2, 0))

    def _us_to_vel(self, vel, acc):
        return vel/7329109, acc/1502
    def _vel_to_us(self, vel, acc):
        return int(7329109*vel), int(1502*acc)
    
    def get_vel_params(self):
        vel, acc = ct.c_int(0), ct.c_int(0)
        self.excmd('ISC_GetVelParams', ct.byref(acc), ct.byref(vel))
        return self._us_to_vel(vel.value, acc.value)
    def set_vel_params(self, vel, acc=10):
        vel, acc = self._vel_to_us(vel, acc)
        self.excmd('ISC_SetVelParams', acc, vel)

    def get_max_positions(self):
        return 360*136533 - 1

    def get_position(self):
        return self.excmd('ISC_GetPosition')
    @rangemethod(0, 'get_max_positions', int)
    def set_position(self, position):
        self.clear_message()
        self.excmd('ISC_MoveToPosition', position)
        self.wait_message(msg=(2, 1))

    def get_angle(self):
        return round(360. * self.get_position() / self.get_max_positions(), 5)
    def set_angle(self, angle):
        position = int(angle/360*self.get_max_positions()) % self.get_max_positions()
        self.set_position(position)

    def set_angle_by(self, angle):
        new_angle = self.get_angle() + angle
        position = int(new_angle/360*self.get_max_positions()) % self.get_max_positions()
        self.set_position(position)