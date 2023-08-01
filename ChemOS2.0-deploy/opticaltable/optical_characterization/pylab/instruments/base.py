import visa as pv
import pyvisa.constants as pv_const
from functools import wraps
from ctypes import create_string_buffer, byref, c_int, c_double, c_uint, Array,  c_long, c_float

## tools
def get_visa_instruments():
    "Get the available visa resources and its corresponding ID."
    visa_lists = pv.ResourceManager().list_resources('?*')
    # visa_lists = pv.ResourceManager('@py').list_resources('?*')
    instrument_list = []
    for v in visa_lists:
        instrument_model = None
        try:
            instrument_v = VisaInstrument(v)
            # instrument_model = instrument_v.test()[:-2]
            instrument_model = instrument_v.manager.get_visa_attribute(3221160169)
            instrument_list.append({
                'visa': v,
                'model': instrument_model.split(',')
            })
        except (pv.VisaIOError):
            instrument_list.append({
                'visa': v,
                'model': ["None"]
            })
    return instrument_list

def list_visa_instruments():
    "Lists the available visa resources and its corresponding instruments."
    instrument_list = get_visa_instruments()

    max_ins = max(instrument_list, key = lambda i: len(i['visa']))
    visa_fmt = '{{:<{:d}s}} >>> {{}}'.format(len(max_ins['visa']))

    for instrument in instrument_list:
        print(visa_fmt.format(
            instrument['visa'],
            ','.join(instrument['model']),
            ))


# ---- helper & decorators ----
class CmdNameMap(object):
    def __init__(self, pairs):
        self.forward = {}
        self.reverse = {}
        for cmd, name in pairs:
            self.forward[cmd] = name.upper()
            self.reverse[name.upper()] = cmd
    def get(self, cmd):
        return self.forward.get(cmd, None)
    def rget(self, name):
        return self.reverse.get(name.upper(), None)
    @property
    def cmds(self):
        return self.forward.keys()
    @property
    def names(self):
        return self.reverse.keys()

def mapsetmethod(cmd_name_map):
    def func_wrapper(func_with_cmd):
        @wraps(func_with_cmd)
        def func_with_name(self, name):
            cmd = cmd_name_map.rget(name)
            if cmd is not None:
                return func_with_cmd(self, cmd)
            else:
                raise ValueError('variable <{}> should be {}'.format(
                    func_with_cmd.__name__[4:], 
                    ', '.join(cmd_name_map.names)))
        return func_with_name
    return func_wrapper

def mapgetmethod(cmd_name_map):
    def func_wrapper(func_ret_cmd):
        @wraps(func_ret_cmd)
        def func_ret_name(self):
            name = cmd_name_map.get(func_ret_cmd(self))
            if name is not None:
                return name
            else:
                raise ValueError('Error in function {}'.format(func_ret_name.__name__))
        return func_ret_name
    return func_wrapper

def rangemethod(min_value, max_value, dtype=None):
    def func_wrapper(func):
# From th260defin.h
        @wraps(func)
        def func_with_range(self, value, min_value=min_value, max_value=max_value):
            # if min, max are strings
            min_value = getattr(self, min_value)() if isinstance(min_value, str) else min_value
            max_value = getattr(self, max_value)() if isinstance(max_value, str) else max_value
            if (dtype is not None) and (not isinstance(value, dtype)):
                raise ValueError('variable <{}> should be of type {}'.format(func.__name__[4:], dtype))
            elif min_value <= value and value <= max_value:
                func(self, value)
            else:
                raise ValueError('variable <{}> should be {} to {}'.format(
                    func.__name__[4:], min_value, max_value ))
        return func_with_range
    return func_wrapper

def add_set_get(Cls):
    def get_method(self, method):
        def err_message(*args):
            print('No method {}'.format(method))
        return getattr(self, method, err_message)

    def set_func(self, **kwargs):
        if len(kwargs.keys()) <= 0:
            methods = dir(self)
            print('set options: {}'.format(', '.join([m[4:] for m in methods if m[:4]=='set_'])))
            return

        # if set(var=value), call set_var(value)
        for var, value in kwargs.items():
            get_method(self, 'set_'+var)(value)

    def get_func(self, *args):
        if len(args) <= 0:
            methods = dir(self)
            print('get options: {}'.format(', '.join([m[4:] for m in methods if m[:4]=='get_'])))
            return

        responses = []
        for var in args:
            value = get_method(self, 'get_'+var)()
            if value is not None:
                responses.append(value)

        if len(responses) == 1:
            return responses[0]
        else:
            return responses

    setattr(Cls, 'set', set_func)
    setattr(Cls, 'get', get_func)
    return Cls

# ctype decorator
class Ref(object):
    def __init__(self, func):
        self.func = func
    def byref(self):
        self.param = self.func()
        return byref(self.param)
    def __getattr__(self, attr):
        if attr == 'value':
            if isinstance(self.param, Array):
                return self.param[:]
            else:
                return self.param.value
def r_str(n):
    return Ref(lambda: create_string_buffer(b"", n))
r_int = lambda: Ref(c_int)
r_uint_arr = lambda n: Ref(c_int * n)
r_double = lambda: Ref(c_double)
r_long = lambda: Ref(c_long)
r_float = lambda: Ref(c_float)

# ctypes wrapper
def ctwrapper(*ctargs):
    def wrap_function(func):
        def ret_func(self, *args):
            refs = []
            wargs = []
            i = 0
            for ctarg in ctargs:
                if isinstance(ctarg, Ref):
                    ref = ctarg
                    refs.append(ref)
                    wargs.append(ref.byref())
                else:
                    cast = ctarg(args[i])
                    wargs.append(cast)
                    i += 1
            func(self, *wargs)

            # decode bytes
            for i, ref in enumerate(refs):
                if isinstance(ref.value, bytes):
                    refs[i] = ref.value.decode('utf-8').rstrip('\x00')
                else:
                    refs[i] = ref.value

            if len(refs) == 0:
                return None
            elif len(refs) == 1:
                return refs[0]
            else:
                return refs
        return ret_func
    return wrap_function

# for wrapping the class
def wrapctypes(list_ctdef):
    def clswrapper(cls):
        for ctdef in list_ctdef:
            def func(self, *args, name = ctdef[0]):
                return getattr(self.lib, name)(*args)
            func = ctwrapper(*ctdef[1:])(func)
            setattr(cls, ctdef[0], func)
        return cls
    return clswrapper





## ---- Base VisaInstruments ----
# rs232 options:
# baud_rate: 9600, 115200
# stop_bits: pv_const.StopBits.one, 
# data_bits: 8,
# parity: pv_const.Parity.none,
# read_termination: '\n', '\r\n'
# write_termination: '\n', '\r\n'
# timeout: <millisecons>

class VisaInstrument(object):
    """
    Base Instrument Manager Class for different instruments
    """

    def __init__(self, visa = None, strip = 0, **kwargs):
        """
        Constructor with optional visa input
        """
        self.manager = None
        self.strip = strip
        if visa:
            # self.set_visa(visa, **kwargs)
            self.manager = pv.ResourceManager().open_resource(visa, **kwargs)
        else:
            self.manager = None

    # def __del__(self):
    #     self.manager.close()

    def test(self):
        "Test the model of manager"
        return str(self.ask('*IDN?'))

    # IO
    def write(self, command):
        "Sends commands to the instrument"
        if command:
            self.manager.write(command)

    def read(self, strip = 0):
        "Reads instrument values"
        if strip > 0:
            return self.manager.read()[:-strip]
        elif self.strip > 0:
            return self.manager.read()[:-self.strip]
        else:
            return self.manager.read()

    def ask(self, query_string):
        "write(query_string) and returns read()"
        if query_string:
            self.write(query_string)
            return self.read()
    def query(self, query_string):
        return self.ask(query_string)






    # ---- set +  get methods ----
    # def set(self, **kwargs):
    #     if len(kwargs.keys()) <= 0:
    #         methods = dir(self)
    #         print('set options: {}'.format(', '.join([m[4:] for m in methods if m[:4]=='set_'])))
    #         return
    #
    #     # if set(var=value), call set_var(value)
    #     for var, value in kwargs.items():
    #         return getattr(self, 'set_'+var)(value)
    #
    # def get(self, *args):
    #     if len(args) <= 0:
    #         methods = dir(self)
    #         print('get options: {}'.format(', '.join([m[4:] for m in methods if m[:4]=='get_'])))
    #         return
    #
    #     for var in args:
    #         return getattr(self, 'get_'+var)()


    # # Methods that must be defined in inherited classes
    # def read_values(self):
    #     "Default read_values() in base instrument class"
    #     return
    # # Load Instrument by Visa
    # def set_visa(self, visa_string, **kwargs):
    #     "Sets the instrument's visa"
    #     if visa_string:






