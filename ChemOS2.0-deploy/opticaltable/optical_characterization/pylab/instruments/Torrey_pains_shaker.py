from . import VisaInstrument, CmdNameMap, mapsetmethod, mapgetmethod, rangemethod, add_set_get
import pyvisa.constants as pv_const
import time


@add_set_get
class TorreyPinesSC20(VisaInstrument):
    def __init__(self, visa, *args, **kwargs):

        # rs232 settings
        self.rs232_settings = {
            'baud_rate': 9600,
            'stop_bits': pv_const.StopBits.one,
            'parity'   : pv_const.Parity.none,
            'data_bits': 8,
            'write_termination': '\r',
            'read_termination': '\r\n',
            'timeout': 5000,
            # 'encoding' : 'utf-8'
        }
        super().__init__(visa, *args, **self.rs232_settings, **kwargs)

    #  # ---- override methods ----
    # def write(self, command):
    #     if command:
    #         self.manager.write(command)
    #         t = self.read()
    #         return t
        self.initialize()

    def initialize(self):
        print(self.query('v'))
        self.get_temparature()

    @rangemethod(-10, 110, int)
    def _set_temparature(self, temparature):
 
        time.sleep(1)
        res = self.query('n%s' %temparature)
        if res == 'e':
            print('SC20:set temparature failed')
        time.sleep(1)
        self.get_temparature_setting()


    def set_temparature(self, temparature, wait = False, tolerance = 1):

        self._set_temparature(temparature)

        if wait:
            print('SC20:Waiting until temparture adjusted')
            while True:
                current_temp = int(self.get_temparature())
                if current_temp >= temparature - tolerance and current_temp <= temparature + tolerance:
                    print('SC20:temparature adjusted to %s' %temparature)
                    break
                time.sleep(10) 


    def get_temparature(self):
        res = self.query('p')
        print('SC20:current temparature:%s' %res)
        return res

    def get_temparature_setting(self):
        res = self.query('s')
        print('SC20:temparature setting:%s' %res)
        return res

    @rangemethod(1, 9, int)
    def mixer_on(self, speed):
        time.sleep(1)
        res = self.query('m%s' %speed)
        if res == 'e':
            print('SC20:set mixer failed')
        time.sleep(1)
        self.get_mixer()

    def mixer_off(self):
        time.sleep(1)
        res = self.query('m%s' %0)
        if res == 'e':
            print('SC20:set mixer off failed')
        time.sleep(1)
        self.get_mixer()

    def get_mixer(self):
        res = self.query('r')
        print('SC20:current mixing speed:%s' %res)
        return res


    
