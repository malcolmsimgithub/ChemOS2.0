from pyvisa.resources import USBInstrument
from pyvisa import ResourceManager
from numpy import frombuffer
from .tlccs_h import ENDPOINT_0_TRANSFERSIZE


# Based on pyvisa USBInstrument
class ModifiedUSBInstrument(object):
    def __init__(self, *args, read_code = 0xC0, write_code = 0x40, **kwargs):
        self.read_code = read_code
        self.write_code = write_code
        kwargs['resource_pyclass'] = USBInstrument

        self.rm = ResourceManager()
        self.manager = self.rm.open_resource(*args, **kwargs)

    def __del__(self):
        self.manager.close()
        self.rm.close()

    def control_out(self, request_type_bitmap_field, request_id, request_value, index, data=""):
        "Modified usb_control_out from the pyvisa repository (just fixing bugs)"
        return self.manager.visalib.usb_control_out(self.manager.session, 
                request_type_bitmap_field, request_id, request_value, index, data)

    ## when future fixes
    # def control_out(self, *args, **kwargs):
    #     return self.manager.control_out(*args, **kwargs)

    def control_in(self, *args, **kwargs):
        return self.manager.control_in(*args, **kwargs)

    def usb_read(self, request_id, request_value, index, length, dtype = None):
        _read_in = lambda v, l: self.control_in(self.read_code, request_id, v, index, l)

        # read in by batches, maximum read length per batch is n_max
        n_max = ENDPOINT_0_TRANSFERSIZE
        data = bytes()
        N, R = length // n_max, length % n_max
        for i in range(N):
            tmp, err = _read_in(request_value+i*n_max, n_max)
            data += tmp
        if R > 0:
            tmp, err = _read_in(request_value+N*n_max, R)
            data += tmp

        # return array if given dtype
        if dtype is not None:
            data = frombuffer(data, dtype=dtype)
        return data, err

    def usb_write(self, request_id, request_value, index, data=""):
        _write_out = lambda v, s, e: self.control_out(self.write_code, request_id, v, index, data[s:e])
        length = len(data)
        if length == 0:
            err = _write_out(request_value, 0, length)
            return err

        # if data length more than 1, then write out by batches
        n_max = ENDPOINT_0_TRANSFERSIZE
        N, R = length // n_max, length % n_max
        for i in range(N):
            err = _write_out(request_value+i*n_max, i*n_max, (i+1)*n_max)
        if R > 0:
            err = _write_out(request_value+N*n_max, N*n_max, N*n_max+R)
        return err

    def set_visa_attribute(self, name, state):
        return self.manager.set_visa_attribute(name, state)

    def get_visa_attribute(self, name):
        return self.manager.get_visa_attribute(name)

    
