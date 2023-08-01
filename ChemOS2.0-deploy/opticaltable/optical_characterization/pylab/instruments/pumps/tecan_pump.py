from .. import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
from .pump import AdvancedPump
import pyvisa.constants as pv_const
import time

MIN_V = 5
MAX_V = 6000
MIN_STOP_V = 50
MAX_STOP_V = 2700
MIN_START_V = 50
MAX_START_V = 1000
MIN_A = 1
MAX_A = 20

CMD_RES_MAP = CmdNameMap([
    (0, 'STANDARD'),
    (1, 'HIGH'), ])


@add_set_get
class XCPump(VisaInstrument, AdvancedPump):

    def __init__(self, visa, addr, init_valve = 12, syringe_volume = 1e-3, ports = {}, info = {}, **kwargs):
        # print(visa, addr, init_valve, syringe_volume, ports)
        self.addr = str(addr)
        self.header = '/' + chr(int(self.addr, 16) + 49)

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'read_termination': '\x03\r\n',
            'timeout': 5000,
        }
        super().__init__(visa, **self.rs232_settings, **kwargs)

        # initialize
        self.init_syringe(init_valve, syringe_volume)
        self.ports = ports
        self.info = info

    def init_syringe(self, init_valve = 12, syringe_volume = 1e-3):

        # set movements
        self.set_resolution('HIGH')
        self.set_acceleration(1) # minimum acceleration
        self.set_start_velocity(50)
        self.set_stop_velocity(50)

        # offset and backlash
        self.write('k0')
        self.write('K96')

        # init
        self.write('Z1,{v},{v}'.format(v=init_valve))

        # syringe volume
        self.syringe_volume = syringe_volume
        self.set_valve_type(12)

    # ---- override methods ----
    def write(self, command):
        return self.ask(command)

    def _ask(self, command):
        self.manager.write(self.header + command + 'R\n')
        return self.read()

    def busy(self):
        _, status = self._ask('Q')
        return (status & 32) == 0

    def ask(self, command):
        if command:
            data, status = self._ask(command)
            while self.busy():
                time.sleep(0.1)
            return data, status

    def read(self):
        response = self.manager.read_raw()
        return response[3:-3], response[2]

    # ---- syringe commands ----
    def get_max_steps(self):
        resolution = self.get_resolution()
        if resolution == 'STANDARD':
            return 3000
        elif resolution == 'HIGH':
            return 24000

    @mapgetmethod(CMD_RES_MAP)
    def get_resolution(self):
        return self.resolution

    @mapsetmethod(CMD_RES_MAP)
    def set_resolution(self, cmd):
        self.write('N{:d}'.format(cmd))
        self.resolution = cmd

    def get_position(self):
        return int(self.ask('?')[0])

    @rangemethod(0, 'get_max_steps')
    def set_position(self, position):
        self.write('A{:.0f}'.format(round(position)))


    # # ---- valve commands ----
    def get_valve_numbers(self):
        num = self.get_valve_type()
        return num

    def get_valve_type(self):
        return self.valve_type

    @rangemethod(1, 12, dtype=int)
    def set_valve_type(self, valve_type):
        self.valve_type = valve_type
        # self.write('U{}')

    def get_valve(self):
        return int(self.ask('?6')[0])

    @rangemethod(1, 'get_valve_numbers', int)
    def set_valve(self, valve):
        self.write('O{}'.format(valve))

    # # ---- motor control commands ----
    def get_velocity(self):
        return int(self.ask('?2')[0])

    @rangemethod(MIN_V, MAX_V)
    def set_velocity(self, velocity):
        self.write('V{:.0f}'.format(velocity))

    def get_start_velocity(self):
        return int(self.ask('?1')[0])

    @rangemethod(MIN_START_V, MAX_START_V)
    def set_start_velocity(self, velocity):
        self.write('v{:.0f}'.format(velocity))

    def get_stop_velocity(self):
        return int(self.ask('?3')[0])

    @rangemethod(MIN_STOP_V, MAX_STOP_V)
    def set_stop_velocity(self, velocity):
        self.write('c{:.0f}'.format(velocity))

    @rangemethod(MIN_A, MAX_A, int)
    def set_acceleration(self, acceleration):
        self.write('L{:d}'.format(acceleration))

