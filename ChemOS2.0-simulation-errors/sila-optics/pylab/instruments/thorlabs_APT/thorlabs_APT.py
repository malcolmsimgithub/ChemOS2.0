from .. import wrapctypes, Ref, r_str, r_int, r_long, r_float, r_double, r_uint_arr
import time, os, sys
import ctypes as ct
from ctypes import byref, c_int, c_long, c_bool, c_float, create_string_buffer
from ctypes.util import find_library

# decode = lambda c: c.value.decode('utf-8')
# strbuf = lambda n: create_string_buffer(b"", n)

list_ctypes_def = [
    ('GetHWSerialNumEx', c_long, c_long, r_long()),
    ('GetNumHWUnitsEx', c_long, r_long()),
    ('GetHWInfo', c_long, r_str(255),c_long, r_str(255), c_long, r_str(255), c_long),
    ('InitHWDevice', c_int),
    ('MOT_EnableHWChannel', c_long),
    ('MOT_DisableHWChannel', c_long),
    ('MOT_GetStatusBits', c_long, r_long()),
    ('MOT_SetHomeParams', c_long, c_long, c_long, c_float, c_float),
    ('MOT_GetHomeParams', c_long, r_long(), r_long(), r_float()),
    ('MOT_SetVelParams', c_long, c_float, c_float, c_float),
    ('MOT_GetVelParams', c_long, r_float(), r_float(), r_float()),
    ('MOT_GetPosition', c_long, r_float()),
    ('MOT_MoveHome', c_long, c_bool),
    ('MOT_MoveAbsoluteEx', c_long, c_float, c_bool),
    ('MOT_MoveRelativeEx', c_long, c_float, c_bool)
]

#dll  = os.path.dirname(os.path.abspath(__file__)

@wrapctypes(list_ctypes_def)
class ThorlabsAPT(object):
    def __init__(self, serial_number):
        filename = "%s/APT.dll" % os.path.dirname(__file__)
        self.lib = ct.windll.LoadLibrary(filename)
        self.serial = serial_number
        print('initializing APT-device')
        self.lib.APTInit()
#        self.lib.EnableEventDlg(False)
        self.lib.APTCleanUp()
        self.lib.APTInit()
        self.InitHWDevice(self.serial)

    def __del__(self):
        self.lib.APTCleanUp()
#        print('closing APT')

    # methods
    def list_available_devices(self):
        """
        Lists all devices connected to the computer.
    
        Returns
        -------
        out : list
            list of available devices. Each device is described by a tuple
            (hardware type, serial number)
        """
        # we have to check for all possible hardware types.
        # Unfortunately I couldn't find a list of all existing hardware types,
        # the list in the C header file is incomplete. Therefore we just check
        # the first 100 type values
        devices = []
        for hwtype in range(100):
            count = self.GetNumHWUnitsEx(hwtype)
                # found an existing hardware type
            if count > 0:
                    # devices are available!!
                    # get their serial number
                for ii in range(count):
                    devices.append((hwtype, self.GetHWSerialNumEx(hwtype, ii)))
        return devices
    
    def hardware_info(self):
        """
        Retrieves hardware information about the devices identified by its
        serial number.
    
        Parameters
        ----------
        serial_number : int
            Serial number identifying the device
    
        Returns
        -------
        out : tuple
            hardware information: (model, software version, hardware notes)
        """
        return self.GetHWInfo(self.serial, 255, 255, 255)

 
    def enable(self):
        """
        Enables the motor (the active channel).
        """
        self.MOT_EnableHWChannel(self.serial)

    def disable(self):
        """
        Disables the motor (the active channel).
        """
        MOT_DisableHWChannel(self.serial)

    def get_velocity_parameters(self):
        """
        Returns current velocity parameters.

        Returns
        -------
        out : tuple
            (minimum velocity, acceleration, maximum velocity)
        """
        return self.MOT_GetVelParams(self.serial)

    def set_velocity_parameters(self, min_vel, accn, max_vel):
        """
        Sets velocity parameters. According to the Thorlabs documentation
        minimum velocity is always 0 and hence is ignored.

        Parameters
        ----------
        min_vel : float
            minimum velocity
        accn : float
            acceleration
        max_vel : float
            maximum velocity
        """
        self.MOT_SetVelParams(self.serial, min_vel, accn, max_vel)

    def get_move_home_parameters(self):
        """
        Returns parameters used when homing.

        Returns
        -------
        out : tuple
            (direction, limiting switch, velocity, zero offset)
        """
        return self.MOT_GetHomeParams(self.serial)

    def set_move_home_parameters(self, direction, lim_switch, velocity, zero_offset):
        """
        Sets parameters used when homing.

        Parameters
        ----------
        direction : int
            home in forward or reverse direction:
            - HOME_FWD = 1 : Home in the forward direction.
            - HOME_REV = 2 : Home in the reverse direction.
        lim_switch : int
            forward limit switch or reverse limit switch:
            - HOMELIMSW_FWD = 4 : Use forward limit switch for home datum.
            - HOMELIMSW_REV = 1 : Use reverse limit switch for home datum.
        velocity : float
            velocity of the motor
        zero_offset : float
            zero offset
        """
        self.MOT_SetHomeParams(self.serial, direction, lim_switch, velocity, zero_offset)

    def move_to(self, value, blocking = False):
        """
        Move to absolute position.

        Parameters
        ----------
        value : float
            absolute position of the motor
        blocking : bool
            wait until moving is finished.
            Default: False
        """
        self.MOT_MoveAbsoluteEx(self.serial, value, blocking)

    def move_by(self, value, blocking = False):
        """
        Move relative to current position.

        Parameters
        ----------
        value : float
            relative distance
        blocking : bool
            wait until moving is finished
            Default: False
        """
        self.MOT_MoveRelativeEx(self.serial, value, blocking)

    @property
    def position(self):
        """
        Position of motor. Setting the position is absolute and non-blocking.
        """
        return self.MOT_GetPosition(self.serial)

    def move_home(self, blocking = False):
        """
        Move to home position.

        Parameters
        ----------
        blocking : bool
            wait until homed
            Default: False
        """
        self.MOT_MoveHome(self.serial, blocking)

#    @property
    def is_in_motion(self):
        """
        Returns whether motor is in motion.
        """
        status_bits = self.MOT_GetStatusBits(self.serial)
        mask = 0x00000010 | 0x00000020 | 0x00000040 | 0x00000080 | 0x00000200
        return bool(status_bits & mask)



