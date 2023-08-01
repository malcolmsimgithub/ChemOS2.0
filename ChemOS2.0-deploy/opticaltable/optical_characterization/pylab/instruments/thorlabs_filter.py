from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time


@add_set_get
class ThorlabsFW212C(VisaInstrument):
    def __init__(self, visa, *args, pos_count=None, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 115200,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r',
            'read_termination': '\r',
            'timeout': 5000,
            # 'encoding' : 'utf-8'
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

        self.initialize(pos_count)

    # def __del__(self):
    #     self.manager.close()

    def device_info(self):
        return self.ask('*idn?')

    def initialize(self, pos_count):
        print(self.device_info())
        self.write('sensors=0')
        self.write('speed=1')
        # print(self.device_info())
        if pos_count:
            self.set_position_count(pos_count)

    # ---- override methods ----
    def write(self, command):
        if command:
            self.manager.write(command)
            t = self.read()
            return t

    # def _ask(self, command):
    #     if command:
    #         print(self.manager.query(command))
    #         return self.read()

    def set_position_count(self, pos_count):
        self.write('pcount={}'.format(pos_count))
    def get_position_count(self):
        return int(self.ask('pcount?'))

    @rangemethod(1, 'get_position_count', int)
    def set_position(self, pos):
        self.write('pos={}'.format(pos))
    def get_position(self):
        return int(self.ask('pos?'))
    
