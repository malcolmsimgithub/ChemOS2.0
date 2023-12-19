from . import VisaInstrument
import pyvisa.constants as pv_const

# TTL-232R-5V-AJ
class TTL232R(VisaInstrument):
    def __init__(self, visa, *args, **kwargs):
        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.even,
            'data_bits': 7,
            'write_termination': '',
            'timeout': 5000,
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

    def on(self):
        self.write(str(0x01))
    def off(self):
        self.write(str(0))

## if there is error with this script, use ftd2xx to set it as bitbang mode and then reuse:
# import ftd2xx as ftd
# d = ftd.open(0)    # Open first FTDI device
# print(d.getDeviceInfo())
# OP = 0x01            # Bit mask for output D0
# d.setBitMode(0x01, 1)  # Set pin as output, and async bitbang mode
