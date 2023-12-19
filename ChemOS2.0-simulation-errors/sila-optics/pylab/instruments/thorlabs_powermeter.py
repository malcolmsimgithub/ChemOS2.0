from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
from pyvisa.resources import USBInstrument
import pyvisa.constants as pv_const

## CONSTANTS##
CMD_MODE_MAP = CmdNameMap([
    ('CURR', 'CURRENT'),
    ('VOLT', 'VOLTAGE'),
    ('POW', 'POWER'), ])

CMD_AUTO_MAP = CmdNameMap([
    ('1', 'ON'),
    ('0', 'OFF'), ])


@add_set_get
class ThorlabsPM100D(VisaInstrument):
    """You may have to change the driver:
        1. Open device manager -> driver -> update driver
        2. Select "USB Test and Measurement Device (IVI)" """

    def __init__(self, visa, averages=None):
        super().__init__(visa, read_termination='\n', resource_pyclass=USBInstrument, timeout=10000)
        if averages:
            self.set_averages(averages)


    # ---- settings ----
    def device_info(self):
        info = self.ask('*IDN?').split(',')
        return '{} {}, serial: {}, firmware: {}'.format(*info)


    @rangemethod(1, 10000, int)
    def set_averages(self, averages):
        self.write('sense:average {}'.format(averages))

    def get_averages(self):
        return int(self.read('sense:average?'))


    @mapsetmethod(CMD_MODE_MAP)
    def set_mode(self, cmd):
        self.write('configure:{}'.format(cmd))

    @mapgetmethod(CMD_MODE_MAP)
    def get_mode(self):
        return self.ask('configure?')


    def measure(self, counts=1):
        value = 0.
        for _ in range(counts):
            value += float(self.ask('read?'))
        return value / counts


    @rangemethod('get_min', 'get_max')
    def set_rang(self, rang):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        self.write('sense:{}:range:upper {}'.format(mode_cmd, rang))

    def get_rang(self):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        return float(self.ask('sense:{}:range:upper?'.format(mode_cmd)))

    def get_max(self):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        return float(self.ask('sense:{}:range:upper? max'.format(mode_cmd)))

    def get_min(self):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        return float(self.ask('sense:{}:range:upper? min'.format(mode_cmd)))

    @mapsetmethod(CMD_AUTO_MAP)
    def set_autorang(self, cmd):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        self.write('sense:{}:range:auto {}'.format(mode_cmd, cmd))

    @mapgetmethod(CMD_AUTO_MAP)
    def get_autorang(self):
        mode_cmd = CMD_MODE_MAP.rget(self.get_mode())
        return self.ask('sense:{}:range:auto?'.format(mode_cmd))


