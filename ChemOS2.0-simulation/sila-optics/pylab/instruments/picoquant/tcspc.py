from .. import wrapctypes, Ref, r_str, r_int, r_double, r_uint_arr
import time, os
import ctypes as ct
from ctypes import byref, c_int, create_string_buffer

# From th260defin.h
LIB_VERSION = "3.1"
MAXDEVNUM = 4
MODE_HIST = 0
MAXLENCODE = 5
MAXINPCHAN = 2
MAXHISTLEN = 32768
FLAG_OVERFLOW = 0x0001

# decode = lambda c: c.value.decode('utf-8')
# strbuf = lambda n: create_string_buffer(b"", n)

list_ctypes_def = [
    ('TH260_GetLibraryVersion', r_str(8)),
    ('TH260_GetHardwareInfo', c_int, r_str(16), r_str(8), r_str(16)),

    ('TH260_OpenDevice', c_int, r_str(8)),
    ('TH260_CloseDevice', c_int),

    ('TH260_Initialize', c_int, c_int),
    ('TH260_GetNumOfInputChannels', c_int, r_int()),

    ('TH260_SetSyncDiv', c_int, c_int),
    ('TH260_SetSyncCFD', c_int, c_int, c_int),
    ('TH260_SetInputCFD', c_int, c_int, c_int, c_int),
    ('TH260_SetSyncChannelOffset', c_int, c_int),
    ('TH260_SetInputChannelOffset', c_int, c_int, c_int),
    ('TH260_SetHistoLen', c_int, c_int, r_int()),
    ('TH260_SetBinning', c_int, c_int),
    ('TH260_SetOffset', c_int, c_int),
    ('TH260_SetStopOverflow', c_int, c_int, c_int),

    ('TH260_GetResolution', c_int, r_double()),
    ('TH260_GetSyncRate', c_int, r_int()),
    ('TH260_GetCountRate', c_int, c_int, r_int()),

    ('TH260_ClearHistMem', c_int),
    ('TH260_StartMeas', c_int, c_int),
    ('TH260_CTCStatus', c_int, r_int()),
    ('TH260_StopMeas', c_int),
    ('TH260_GetHistogram', c_int, r_uint_arr(MAXHISTLEN), c_int, c_int),
]

module_path = os.path.dirname(os.path.abspath(__file__))

@wrapctypes(list_ctypes_def)
class TH260(object):
    def __init__(self, device = 0, dll=os.path.join(module_path, "th260lib64.dll")):
        self.lib = ct.CDLL(dll)
        self.dev = device

        self.serial = self.TH260_OpenDevice(self.dev)
        self.TH260_Initialize(self.dev, MODE_HIST)

        self.set_histlen(MAXLENCODE)

    # methods
    def close(self, i):
        self.TH260_CloseDevice(i)
    def close_all(self):
        for i in range(0, MAXDEVNUM):
            self.close(i)

    def get_lib_version(self):
        return self.TH260_GetLibraryVersion()
    def get_hardware_info(self):
        return self.TH260_GetHardwareInfo(self.dev)
    def get_numchannels(self):
        return self.TH260_GetNumOfInputChannels(self.dev)

    # settings
    def set_sync_divider(self, sync_divider):
        self.TH260_SetSyncDiv(self.dev, sync_divider)

    def set_cfd(self, channel, cfd_level, zero_cross):
        if channel == 'SYNC':
            self.TH260_SetSyncCFD(self.dev, cfd_level, zero_cross)
        else:
            self.TH260_SetInputCFD(self.dev, channel, cfd_level, zero_cross)

    def set_channel_offset(self, channel, offset):
        if channel == 'SYNC':
            self.TH260_SetSyncChannelOffset(self.dev, offset)
        else:
            self.TH260_SetInputChannelOffset(self.dev, channel, offset)

    def set_histlen(self, histlen):
        return self.TH260_SetHistoLen(self.dev, histlen)
    def set_binning(self, binning): # need to update this function (convert binning to something meaningful)
        self.TH260_SetBinning(self.dev, binning)
    def set_offset(self, offset):
        self.TH260_SetOffset(self.dev, offset)
   # def set_overflow(self, flag, overflow): # need to change this function (flag=0 means no stop
   #     self.TH260_SetStopOverflow(self.dev, flag, overflow) 
    def set_overflow(self, overflow):         
        if overflow == 0:
            self.TH260_SetStopOverflow(self.dev, 0, overflow)
        else:
            self.TH260_SetStopOverflow(self.dev, 1, overflow)


    def get_resolution(self):
        return self.TH260_GetResolution(self.dev)
    def get_sync_rate(self):
        return self.TH260_GetSyncRate(self.dev)
    def get_count_rate(self, channel):
        return self.TH260_GetCountRate(self.dev, channel)

    # meausrement
    def measure(self, total_time, channel=0):

        self.TH260_ClearHistMem(self.dev)
        self.TH260_StartMeas(self.dev, int(total_time*1e3))
        while self.TH260_CTCStatus(self.dev) == 0:
            time.sleep(0.1)
        self.TH260_StopMeas(self.dev) 
        return self.TH260_GetHistogram(self.dev, channel, 0).copy()


if __name__ == "__main__":
    pass



# obsolete (for wrapper based)
    # TH260 methods
    # def __getattr__(self, func_name):
    #     if func_name[:6] == 'TH260_':
    #         func = getattr(self.lib, func_name, None)
    #         if func is None:
    #             raise ValueError('{:s} not found'.format(func_name))
    #         # create func
    #         def func_errcheck(*args):
    #             retcode = func(*args)
    #             message = strbuf(40)
    #             if retcode < 0:
    #                 self.lib.TH260_GetErrorString(message, c_int(retcode))
    #                 self.close_all()
    #                 raise ValueError("{:s} error {:d} ({:s})".format(func_name, retcode, decode(message)))
    #             # return retcode, decode(message)
    #         return func_errcheck
    #     else:
    #         raise ValueError('{} not found'.format(func_name))

# class Ref(object):
#     def __init__(self, func):
#         self.func = func
#     def __call__(self):
#         return self.func()
#
# def refbuf(n):
#     return Ref(lambda: strbuf(n))
#
# def ctwrapper(*ctargs):
#     def wrap_function(func):
#         def ret_func(*args):
#             refs = []

#             wargs = []
#             i = 0
#             for ctarg in ctargs:
#                 if isinstance(ctarg, Ref):
#                     ref = ctarg()
#                     refs.append(ref)
#                     wargs.append(ref)
#                 else:
#                     cast = ctarg(args[i])
#                     wargs.append(cast)
#                     i += 1
#             func(*wargs)
#
#             # decode bytes
#             for i, ref in enumerate(refs):
#                 if isinstance(ref.value, bytes):
#                     refs[i] = ref.value.decode('utf8')
#
#             if len(refs) == 0:
#                 return None
#             elif len(refs) == 1:
#                 return refs[0]
#             else:
#                 return refs
#         return ret_func
#     return wrap_function
if __name__ == "__main__":
    pass
