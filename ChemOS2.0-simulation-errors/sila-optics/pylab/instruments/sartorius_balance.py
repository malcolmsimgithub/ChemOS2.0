from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time


@add_set_get
class EntrisBalance(VisaInstrument):
    def __init__(self, visa, *args, **kwargs):
        self.rs232 = {
            'baud_rate': 19200,
            'parity'   : pv_const.Parity.odd,
            'stop_bits': pv_const.StopBits.one,
            'data_bits': 7,
            'read_termination': '\r\n',
            'write_termination': '\r\n',
            'timeout': 15000,
        }
        super().__init__(visa, *args, **self.rs232, **kwargs)

    def __del__(self):
        print('Balance closed')
        self.manager.close()

    def device_info(self):
        return 'Sartoris balance model: {}'.format(self.ask('\x1bx1_').replace(' ', ''))

    def tare(self):
        # may need to wait
        self.write('\x1bT')

    def get_weight(self):
        print_str = self.ask('\x1bP')
        return float(print_str[6:16].replace(' ', ''))

    
