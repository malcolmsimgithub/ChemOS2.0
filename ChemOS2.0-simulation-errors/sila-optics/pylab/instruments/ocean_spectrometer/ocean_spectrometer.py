# import clr
import os, sys
import numpy as np
import time
from typing import *
# from pylab.instruments.ocean_spectrometer import seabreeze

# from pylab.instruments.ocean_spectrometer.seabreeze.spectrometers import Spectrometer, list_devices
# from pylab.instruments.ocean_spectrometer.seabreeze.cseabreeze import SeaBreezeAPI
from seabreeze.spectrometers import Spectrometer, list_devices
from seabreeze.cseabreeze import SeaBreezeAPI


# device = list_devices()


def list_device():
    for device in list_devices():
        print(device)

class OceanSpectrometer():

    def __init__(self, serial_number = None, verbose = False):

        self.verbose = verbose

        if serial_number:
            self.spec = Spectrometer.from_serial_number(serial_number)
        else:
            self.spec = Spectrometer.from_first_available()

        self.wavelengths = self.spec.wavelengths
        self.minimum_integration_time = self.spec.minimum_integration_time_micros/1e6
        self.maxinum_intensity = self.spec.max_intensity
        
        print(self.spec)


    def set_integration_time(self, integration_time: int):
        """
        set integration_time

        Parameters
        ------------------------------------
        integration_time : int
            integration time in seconds
        """
        if integration_time < self.minimum_integration_time:
            raise Exception('integration time should be more than %s s' %self.minimum_integration_time)

        self.spec.integration_time_micros(integration_time*1E6)
        time.sleep(0.5)
        if self.verbose:
            print('%s: integration time is set to %s s' %(self.spec.model, integration_time))


    def measure(self, repeats: int = 1, integration_time: int = None, dark_correction: bool = True, nonlinearity_correction: bool = True, normalize: bool = True):

        if integration_time is not None:
            self.set_integration_time(integration_time)

        data = []
        saturation = False
        self.spec.f.data_buffer.clear()
        for i in range(repeats):
            data.append(self.spec.intensities(correct_dark_counts = dark_correction, correct_nonlinearity= nonlinearity_correction))
            if max(data[i]) > self.maxinum_intensity:
                saturation = True
            if normalize:
                data[i] /= self.maxinum_intensity
        
        if saturation:
            print('%s: Caution! Spectrometer is saturated.')
        
        return np.average(np.asarray(data), axis = 0)
        

    def measure_trace(self, repeats = 1, integration_time = None, dark_correction = True, nonlinearity_correction = True, normalize = True):
        
        if integration_time is not None:
            self.set_integration_time(integration_time)

        time_list = []
        data = []
        saturation = False

        self.spec.f.data_buffer.clear()
        for i in range(repeats):
            time_list.append(time.time())
            if i == 0:
                time_start = time_list[0]
            time_list[i] -= time_start
            data.append(self.spec.intensities(correct_dark_counts = dark_correction, correct_nonlinearity= nonlinearity_correction))
            if max(data[i]) > self.maxinum_intensity:
                saturation = True
            if normalize:
                data[i] /= self.maxinum_intensity
                
        if saturation:
            print('%s: Caution! Spectrometer is saturated.')

        results = {'time' : time_list, 'average' : np.average(np.asarray(data), axis = 0), 'time_trace' : np.asarray(data)}

        return results




# API = SeaBreezeAPI()

# Device = API.list_devices()


# spec = Device[0]

# print(spec.serial_number)
# print(spec.features)

# spec.f.spectrometer.set_integration_time_micros(int(1000000))

# spec = Spectrometer.from_serial_number('QEP03644')

# spec.integration_time_micros(100000)
# time.sleep(1)
# print(spec.features)

# wavelengths = spec.wavelengths()

# spec1 = spec.intensities()

# time.sleep(10)

# spec2 = spec.intensities()

# t = spec1/spec2

# plt.plot(wavelengths,spec1)
# plt.plot(wavelengths,spec2)
# plt.show()

# plt.plot(wavelengths,t)
# plt.show()

# spec = Spectrometer.from_first_available()
# print(spec.serial_number)
# from ctypes import windll

# #sys.path.append('C:/Program Files(x86)/Chromeleon/SDK Documentation/Assemblies')
# sys.path.append('C:\Program Files\Ocean Optics\OmniDriverSPAM\OOI_HOME')
# # sys.path.append('C:\Program Files\Ocean Optics\OmniDriverSPAM\_jvm\bin')
# # clr.AddReference('NETOmniDriver-NET40')
# # clr.AddReference('System')

# # from System import Uri, String, Transactions, TimeSpan, Double, Nullable, Int32

# print(os.path.exists('C:\Program Files\Ocean Optics\OmniDriverSPAM\OOI_HOME\OmniDriver32.dll'))

# dll = windll.LoadLibrary('C:/Program Files/Ocean Optics/OmniDriverSPAM/OOI_HOME/OmniDriver64.dll')
# from OmniDriver import NETWrapper
# import OmniDriver

# wrapper = NETWrapper()
# a, b = wrapper.getName(Int32(0), Int32(0))

# wrapper.GetType()

# print(dir(wrapper))

# wrapper.closeAllSpectrometers()

# wrapper = NETWrapper


# # print(dir(wrapper))

# # wrapper.getDetector()
# num = wrapper.openAllSpectrometers()



# from Thermo.Chromeleon.Sdk.Common import CmSdkScope, CmSdk
# from Thermo.Chromeleon.Sdk.Interfaces import Data #, UserManagement

# from System import Uri, String, Transactions, TimeSpan, Double, Nullable

# class OceanSpectrometer():

#     def __init__(self):
        
#         self.MSReader = Thermo_MSReader()
#         self.Analyzer = Chrom_Analyzer()
#         self.Device = None
#         if not offline:
#             from .ThermoFisher_Chromeleon import Chromeleon
#             if Chemspeed:
#                 self.Device = Chromeleon("chrom://fnsvnr2/ChromeleonLocal/Instrument Data", 'UofT_LCMS_QE_Chemspeed')
#             else:
#                 self.Device = Chromeleon("chrom://fnsvnr2/ChromeleonLocal/Instrument Data", 'UofT_LCMS_QE')
#         if DB:
#             from .MDB_client import MDB_client_HPLCMS
#             self.db = MDB_client_HPLCMS(config_file = '%s\.config_MDB_client.dat' %filedir, admin=True)
#         self.target_dict_file = '%s/MS_target_dict' %filedir
#         self.target_dict = self.load_target_dict()
#         #self.target_type = ['BA', 'BMIDA', 'X-BMIDA', 'product', 'others']
#         self.target_type = ['BA', 'BMIDA', 'X-BMIDA', 'product', 'product_radical','others']


#     def instrument_method(func):
#         def wrapper(self, *args, **kwargs):
#             if not self.Device:
#                 raise Exception ('Function <{}> needs instruments'.format(func.__name__))
#             else: 
#                 func(self, *args, **kwargs)
#         return wrapper

