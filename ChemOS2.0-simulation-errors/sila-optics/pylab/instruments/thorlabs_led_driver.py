from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const



NUM_LEDS = 4
LED_BITS = [0x20, 0x80, 0x200, 0x800]

CMD_SEL_MAP = CmdNameMap([
    ('0', 'MULTI'),
    ('1', 'SINGLE'), ])

CMD_MODE_MAP = CmdNameMap([
    ('0', 'CURRENT'),
    ('1', 'PERCENTAGE'),
    ('2', 'EXTERNAL'), ])

CMD_STATE_MAP = CmdNameMap([
    ('0', 'OFF'),
    ('1', 'ON'), ])


@add_set_get
class ThorlabsLED(object):
    "Don't create this directly, call from DC4100"
    
    def __init__(self, channel, manager):
        if 0 <= channel and channel < NUM_LEDS:
            self.channel = channel
            self.parent = manager
        else:
            raise ValueError('<channel> should be between 0~{}'.format(NUM_LEDS-1))
        self.set_limit(self.get_max_limit())

    def device_info(self):
        return 'Thorlabs LED {}, serial: {}, wavelength: {} nm'.format(
            self.ask('hn? {c}'), 
            self.ask('hs? {c}'), 
            self.get_wavelength())

    # ---- override ----
    # format the command with channel where {c}
    def write(self, command):
        self.parent.write(command.format(c=self.channel))
    def read(self): 
        return self.parent.read()
    def ask(self, command): 
        return self.parent.ask(command.format(c=self.channel))

    # ---- led operations ----
    @rangemethod(0, 'get_limit')
    def set_current(self, current):
        "units should be in A"
        self.write('cc {{c}} {:d}'.format(round(1e3 * current)))
        # if current > 0:
        #     self.on()
        # else:
        #     self.off()

    def get_current(self):
        return 1e-3 * float(self.ask('cc? {c}'))

    @rangemethod(0, 'get_max_limit')
    def set_limit(self, limit):
        "units should be in A"
        self.write('l {{c}} {:d}'.format(round(1e3 * limit)))
        
    def get_limit(self):
        return 1e-3 * float(self.ask('l? {c}'))

    def get_max_limit(self):
        return 1e-3 * float(self.ask('ml? {c}'))

    @rangemethod(0, 100)
    def set_percentage(self, percentage):
        self.write('bp {{c}} {}'.format(percentage))
        # if percentage > 0:
        #     self.on()
        # else:
        #     self.off()

    def get_percentage(self):
        return float(self.ask('bp? {c}'))

    def get_wavelength(self):
        return float(self.ask('wl? {c}'))

    @mapsetmethod(CMD_STATE_MAP)
    def set_state(self, cmd):
        self.write('o {{c}} {}'.format(cmd))

    @mapgetmethod(CMD_STATE_MAP)
    def get_state(self):
        return self.ask('o? {c}')

    def on(self):
        self.set_state('ON')
    def off(self):
        self.set_state('OFF')








@add_set_get
class ThorlabsDC4100(VisaInstrument):
    def __init__(self, visa, *args, **kwargs):
        rs232_settings = {
            'baud_rate': 115200,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'read_termination': '\r\n',
        }
        # self.manager = VisaInstrument(visa, *args, **rs232_settings, **kwargs)
        super().__init__(visa, *args, **rs232_settings, **kwargs)
        self.set_selection('MULTI')

        # load led
        self.led = {}
        self.load_led()

    def __del__(self):
        # for _, led in self.led.items():
        #     del led
        self.manager.close()

    def load_led(self):
        status = self.get_status()
        for i in range(NUM_LEDS):
            if (status & LED_BITS[i]) == 0:
                self.led[i+1] = ThorlabsLED(i, self)

    # ---- settings ----
    def device_info(self):
        return '{} {}, serial: {}, firmware: {}'.format(
            self.ask('h?'), 
            self.ask('n?'),
            self.ask('s?'),
            self.ask('v?'))

    def get_status(self):
        return int(self.ask('r?'))

    @mapsetmethod(CMD_SEL_MAP)
    def set_selection(self, cmd):
        self.write('sm {}'.format(cmd))

    @mapgetmethod(CMD_SEL_MAP)
    def get_selection(self):
        return self.ask('sm?')

    @mapsetmethod(CMD_MODE_MAP)
    def set_mode(self, cmd):
        self.write('o -1 0')
        self.write('m {}'.format(cmd))

    @mapgetmethod(CMD_MODE_MAP)
    def get_mode(self):
        return self.ask('m?')


    # def set_mode(self, mode):
    #     cmd = MODE_TO_CMD.get(mode.upper(), None)
    #     if cmd is not None:
    #         self.write('o -1 0')
    #         self.write('m {}'.format(cmd))
    #     else:
    #         raise ValueError('<mode> should be {}'.format(', '.join(MODE_TO_CMD.keys())))
    #
    # def get_mode(self):
    #     mode = CMD_TO_MODE.get(self.ask('m?'), None)
    #     if mode is not None:
    #         return mode

    # def set_selection(self, selection):
    #     cmd = SEL_TO_CMD.get(selection.upper(), None)
    #     if cmd is not None:
    #         self.write('sm {}'.format(cmd))
    #     else:
    #         raise ValueError('<selection> should be {}'.format(SEL_TO_CMD.keys()))
    #
    # def get_selection(self):
    #     selection = CMD_TO_SEL.get(self.ask('sm?'), None)
    #     if selection is not None:
    #         return selection

    # def set_state(self, state):
    #     cmd = STATE_TO_CMD.get(state.upper(), None)
    #     if cmd is not None:
    #         self.write('o {{c}} {}'.format(cmd))
    #     else:
    #         raise ValueError('<state> should be {}'.format(', '.join(STATE_TO_CMD.keys())))
    #
    # def get_state(self):
    #     state = CMD_TO_STATE.get(self.ask('o? {c}'), None)
    #     if state is not None:
    #         return state



    # def set_percentage(self, percentage):
    #     if 0<=percentage and percentage <= 100:
    #         self.write('bp {{c}} {}'.format(percentage))
    #     else:
    #         raise ValueError('percentage should be 0~100')
    #
    #     if percentage > 0:
    #         self.on()
    #     else:
    #         self.off()

    # def set_limit(self, limit):
    #     "Limit should be in A"
    #     max_limit = self.get_max_limit()
    #     if 0<=limit and limit <= max_limit:
    #         self.write('l {{c}} {:d}'.format(round(1e3 * limit)))
    #     else:
    #         raise ValueError('limit should be 0~{} A'.format(max_limit))
    # def set_current(self, current):
    #     "Current should be in A"
    #     limit = self.get_limit()
    #     if 0<=current and current <= limit:
    #         self.write('cc {{c}} {:d}'.format(round(1e3 * current)))
    #     else:
    #         raise ValueError('current should be 0~{} A'.format(limit))
    #
    #     if current > 0:
    #         self.on()
    #     else:
    #         self.off()
