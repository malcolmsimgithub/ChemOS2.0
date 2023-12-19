#### using Kinesis
from .. import CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
from ctypes import cdll, create_string_buffer
import ctypes as ct
from ctypes.wintypes import DWORD
import os, time, sys # clr

__all__ = ['ThorlabsKSC101']


CMD_MODE_MAP = CmdNameMap([
    (0x01, 'MANUAL'),
    (0x02, 'SINGLE'),
    (0x03, 'AUTO'),
    (0x04, 'TRIGGERED'), 
])

CMD_STATE_MAP = CmdNameMap([
    (0x01, 'ACTIVE'),
    (0x02, 'INACTIVE'),
])

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
class ThorlabsKSC101(object):
    def __init__(self, visa, dllpath=dll_path, 
            dllname=r'Thorlabs.MotionControl.KCube.Solenoid.dll'):

        os.environ['PATH'] = dllpath + os.pathsep + os.environ['PATH']
        self.apt = cdll.LoadLibrary(dllname)

        # check visa num
        if self.apt.TLI_BuildDeviceList() != 0:
            raise ValueError('Something wrong with device')
        self.apt.TLI_GetDeviceListSize()
        data = create_string_buffer(100)
        self.apt.TLI_GetDeviceListByTypeExt(data, 100, 68)
        if visa in str(data.value):
            self.serial = ct.c_char_p(str.encode(visa))
            self.excmd('SC_Open')
        else:
            self.serial = ct.c_char_p(str.encode(''))
            print('Device Not Found')

    def __del__(self):
        self.excmd('SC_Close')

    def excmd(self, command, *args):
        resp = getattr(self.apt, command)(self.serial, *args)
        time.sleep(0.02)
        return resp

    def device_info(self):
        info = DeviceInfo()
        self.apt.TLI_GetDeviceInfo(self.serial, ct.byref(info))
        return '{}, serial: {}'.format(info.description.decode(), info.serialNo.decode())


    @mapsetmethod(CMD_MODE_MAP)
    def set_mode(self, cmd):
        self.excmd('SC_SetOperatingMode', cmd)

    @mapgetmethod(CMD_MODE_MAP)
    def get_mode(self):
        return self.excmd('SC_GetOperatingMode')

    @mapsetmethod(CMD_STATE_MAP)
    def set_state(self, cmd):
        self.excmd('SC_SetOperatingState', cmd)

    @mapgetmethod(CMD_STATE_MAP)
    def get_state(self):
        return self.excmd('SC_GetOperatingState')

    def open(self):
        self.set_state('ACTIVE')
    def close(self):
        self.set_state('INACTIVE')


