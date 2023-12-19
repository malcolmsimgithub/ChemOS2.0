from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time

MIN_T = 20
MAX_T = 300
MIN_R = 100
MAX_R = 1400

CMD_STATE_MAP = CmdNameMap([
    ('START', 'ON'),
    ('STOP', 'OFF'), ])

@add_set_get
class HotplateHei(VisaInstrument):
    def __init__(self, visa, *args, **kwargs):
        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.even,
            'data_bits': 7,
            'read_termination': '\r\n',
            'timeout': 5000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        # init
        # self.reset()
        self.write('PA_NEW')

    def device_info(self):
        return 'Heidolph, firmware: {}'.format(self.ask('SW_VERS')[9:])

    def reset(self):
        self.write('RESET')

    # ---- override methods ----
    def write(self, command):
        if command:
            time.sleep(0.1)
            self.manager.write(command)
            return self.read()

    def ask(self, command):
        return self.write(command)

    
    # ---- heating ----
    @rangemethod(MIN_T, MAX_T)
    def set_temperature(self, temperature):
        self.write('OUT_SP_1 {}'.format(temperature))
        if temperature > MIN_T:
            self.t_on()
        else:
            self.t_off()

    def get_temperature(self):
        return float(self.ask('IN_PV_4')[8:])

    def set_T(self, temperature):
        self.set_temperature(temperature)
    def get_T(self):
        return self.get_temperature()

    @mapsetmethod(CMD_STATE_MAP)
    def set_temperature_state(self, cmd):
        self.write('{}_1'.format(cmd))

    # @mapgetmethod(CMD_STATE_MAP)
    # def get_state(self):
    #     return self.ask('o? {c}')

    def t_on(self):
        self.set_temperature_state('ON')
    def t_off(self):
        self.set_temperature_state('OFF')

    # ---- spinning ----
    @rangemethod(MIN_R, MAX_R)
    def set_rotation(self, rotation):
        self.write('OUT_SP_3 {}'.format(rotation))
        if rotation > MIN_R:
            self.r_on()
        else:
            self.r_off()

    def get_rotation(self):
        return float(self.ask('IN_PV_5')[8:])

    def set_R(self, rotation):
        self.set_rotation(rotation)
    def get_R(self):
        return self.get_rotation()

    @mapsetmethod(CMD_STATE_MAP)
    def set_rotation_state(self, cmd):
        self.write('{}_2'.format(cmd))

    def spin(self, rotation):
        self.set_rotation(rotation)
        while self.get_rotation() < rotation*0.95:
            time.sleep(0.2)

    # @mapgetmethod(CMD_STATE_MAP)
    # def get_state(self):
    #     return self.ask('o? {c}')

    def r_on(self):
        self.set_rotation_state('ON')
    def r_off(self):
        self.set_rotation_state('OFF')


# class HotplateHei(VisaInstrument):
#     def __init__(self, visa, *args, **kwargs):
#         # rs232 settings
#         self.rs232_settings = {
#             'baud_rate': 9600,
#             'stop_bits': pv_const.StopBits.one,
#             'parity'   : pv_const.Parity.even,
#             'data_bits': 7,
#             'read_termination': '\r\n',
#             'timeout': 5000,
#         }
#         super().__init__(visa, *args, **self.rs232_settings, **kwargs)
#         # init
#         self.write('PA_NEW')
#
#     def device_info(self):
#         return 'Heidolph, firmware: {}'.format(self.ask('SW_VERS')[9:])
#
#     def reset(self):
#         self.write('RESET')
#
#     # ---- override methods ----
#     def write(self, command):
#         if command:
#             time.sleep(0.1)
#             self.manager.write(command)
#             return self.read()
#
#     def ask(self, command):
#         return self.write(command)
#
#     
#     # ---- heating ----
#     @rangemethod(MIN_T, MAX_T)
#     def set_temperature(self, temperature):
#         self.write('OUT_SP_1 {}'.format(temperature))
#         if temperature > MIN_T:
#             self.write('START_1')
#         elif temperature < MIN_T:
#             self.write('OUT_SP_1 {}'.format(MIN_T))
#             self.write('STOP_1')
#         else:
#             raise ValueError('<temperature> should be smaller than {:d}'.format(MAX_T))
#
#     def get_temperature(self):
#         return float(self.ask('IN_PV_4')[8:])
#
#     # ---- spinning ----
#     def set_rotation(self, rotation):
#         if MIN_R <= rotation and rotation <= MAX_R:
#             print(rotation)
#             self.write('OUT_SP_3 {}'.format(rotation))
#             self.write('START_2')
#         elif rotation < MIN_R:
#             self.write('OUT_SP_3 {}'.format(MIN_R))
#             self.write('STOP_2')
#         else:
#             raise ValueError('<rotation> should be smaller than {:d}'.format(MAX_R))
#
#     def get_rotation(self):
#         return float(self.ask('IN_PV_5')[8:])
#
#
#
#
#
