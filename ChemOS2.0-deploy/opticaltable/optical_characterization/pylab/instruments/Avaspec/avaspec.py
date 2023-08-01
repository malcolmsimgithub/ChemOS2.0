import sys
import ctypes
import ctypes.wintypes
import struct
import globals
from PyQt5.QtCore import *

AVS_SERIAL_LEN = 10
VERSION_LEN = 16
USER_ID_LEN = 64
WM_MEAS_READY = 0x8001

dll = "C:/Users/Matterlab/Dropbox/PythonLab/pylab/instruments/Avaspec/avaspecx64.dll"

class AvsIdentityType(ctypes.Structure):
  _pack_ = 1
  _fields_ = [("SerialNumber", ctypes.c_char * AVS_SERIAL_LEN),
              ("UserFriendlyName", ctypes.c_char * USER_ID_LEN),
              ("Status", ctypes.c_char)]

class MeasConfigType(ctypes.Structure):
  _pack_ = 1
  _fields_ = [("m_StartPixel", ctypes.c_uint16),
              ("m_StopPixel", ctypes.c_uint16),
              ("m_IntegrationTime", ctypes.c_float),
              ("m_IntegrationDelay", ctypes.c_uint32),
              ("m_NrAverages", ctypes.c_uint32),
              ("m_CorDynDark_m_Enable", ctypes.c_uint8), # nesting of types does NOT work!!
              ("m_CorDynDark_m_ForgetPercentage", ctypes.c_uint8),
              ("m_Smoothing_m_SmoothPix", ctypes.c_uint16),
              ("m_Smoothing_m_SmoothModel", ctypes.c_uint8),
              ("m_SaturationDetection", ctypes.c_uint8),
              ("m_Trigger_m_Mode", ctypes.c_uint8),
              ("m_Trigger_m_Source", ctypes.c_uint8),
              ("m_Trigger_m_SourceType", ctypes.c_uint8),
              ("m_Control_m_StrobeControl", ctypes.c_uint16),
              ("m_Control_m_LaserDelay", ctypes.c_uint32),
              ("m_Control_m_LaserWidth", ctypes.c_uint32),
              ("m_Control_m_LaserWaveLength", ctypes.c_float),
              ("m_Control_m_StoreToRam", ctypes.c_uint16)]

class DeviceConfigType(ctypes.Structure):
  _pack_ = 1
  _fields_ = [("m_Len", ctypes.c_uint16),
              ("m_ConfigVersion", ctypes.c_uint16),
              ("m_aUserFriendlyId", ctypes.c_char * USER_ID_LEN),
              ("m_Detector_m_SensorType", ctypes.c_uint8),                      
              ("m_Detector_m_NrPixels", ctypes.c_uint16),
              ("m_Detector_m_aFit", ctypes.c_float * 5),
              ("m_Detector_m_NLEnable", ctypes.c_bool),
              ("m_Detector_m_aNLCorrect", ctypes.c_double * 8),
              ("m_Detector_m_aLowNLCounts", ctypes.c_double),
              ("m_Detector_m_aHighNLCounts", ctypes.c_double),
              ("m_Detector_m_Gain", ctypes.c_float * 2),
              ("m_Detector_m_Reserved", ctypes.c_float),
              ("m_Detector_m_Offset", ctypes.c_float * 2),
              ("m_Detector_m_ExtOffset", ctypes.c_float),
              ("m_Detector_m_DefectivePixels", ctypes.c_uint16 * 30),
              ("m_Irradiance_m_IntensityCalib_m_Smoothing_m_SmoothPix", ctypes.c_uint16),
              ("m_Irradiance_m_IntensityCalib_m_Smoothing_m_SmoothModel", ctypes.c_uint8),
              ("m_Irradiance_m_IntensityCalib_m_CalInttime", ctypes.c_float),
              ("m_Irradiance_m_IntensityCalib_m_aCalibConvers", ctypes.c_float * 4096),
              ("m_Irradiance_m_CalibrationType", ctypes.c_uint8),
              ("m_Irradiance_m_FiberDiameter", ctypes.c_uint32),  
              ("m_Reflectance_m_Smoothing_m_SmoothPix", ctypes.c_uint16),
              ("m_Reflectance_m_Smoothing_m_SmoothModel", ctypes.c_uint8),
              ("m_Reflectance_m_CalInttime", ctypes.c_float),
              ("m_Reflectance_m_aCalibConvers", ctypes.c_float * 4096),
              ("m_SpectrumCorrect", ctypes.c_float * 4096),
              ("m_StandAlone_m_Enable", ctypes.c_bool),
              ("m_StandAlone_m_Meas_m_StartPixel", ctypes.c_uint16),
              ("m_StandAlone_m_Meas_m_StopPixel", ctypes.c_uint16),
              ("m_StandAlone_m_Meas_m_IntegrationTime", ctypes.c_float),
              ("m_StandAlone_m_Meas_m_IntegrationDelay", ctypes.c_uint32),
              ("m_StandAlone_m_Meas_m_NrAverages", ctypes.c_uint32),
              ("m_StandAlone_m_Meas_m_CorDynDark_m_Enable", ctypes.c_uint8), 
              ("m_StandAlone_m_Meas_m_CorDynDark_m_ForgetPercentage", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_Smoothing_m_SmoothPix", ctypes.c_uint16),
              ("m_StandAlone_m_Meas_m_Smoothing_m_SmoothModel", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_SaturationDetection", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_Trigger_m_Mode", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_Trigger_m_Source", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_Trigger_m_SourceType", ctypes.c_uint8),
              ("m_StandAlone_m_Meas_m_Control_m_StrobeControl", ctypes.c_uint16),
              ("m_StandAlone_m_Meas_m_Control_m_LaserDelay", ctypes.c_uint32),
              ("m_StandAlone_m_Meas_m_Control_m_LaserWidth", ctypes.c_uint32),
              ("m_StandAlone_m_Meas_m_Control_m_LaserWaveLength", ctypes.c_float),
              ("m_StandAlone_m_Meas_m_Control_m_StoreToRam", ctypes.c_uint16),
              ("m_StandAlone_m_Nmsr", ctypes.c_int16),
              ("m_StandAlone_m_Reserved", ctypes.c_uint8 * 12), # SD Card, do not use
              ("m_Temperature_1_m_aFit", ctypes.c_float * 5),
              ("m_Temperature_2_m_aFit", ctypes.c_float * 5),
              ("m_Temperature_3_m_aFit", ctypes.c_float * 5),
              ("m_TecControl_m_Enable", ctypes.c_bool),
              ("m_TecControl_m_Setpoint", ctypes.c_float),
              ("m_TecControl_m_aFit", ctypes.c_float * 2),
              ("m_ProcessControl_m_AnalogLow", ctypes.c_float * 2),
              ("m_ProcessControl_m_AnalogHigh", ctypes.c_float * 2),
              ("m_ProcessControl_m_DigitalLow", ctypes.c_float * 10),
              ("m_ProcessControl_m_DigitalHigh", ctypes.c_float * 10),
              ("m_EthernetSettings_m_IpAddr", ctypes.c_uint32),
              ("m_EthernetSettings_m_NetMask", ctypes.c_uint32),
              ("m_EthernetSettings_m_Gateway", ctypes.c_uint32), 
              ("m_EthernetSettings_m_DhcpEnabled", ctypes.c_uint8), 
              ("m_EthernetSettings_m_TcpPort", ctypes.c_uint16),
              ("m_EthernetSettings_m_LinkStatus", ctypes.c_uint8),  
              ("m_Reserved", ctypes.c_uint8 * 9720),
              ("m_OemData", ctypes.c_uint8 * 4096)]

def AVS_Init(x):
    import os   
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int)
    paramflags = (1, "port",),
    AVS_Init = prototype(("AVS_Init", lib), paramflags)
    ret = AVS_Init(x) 
    return ret 

def AVS_GetNrOfDevices():
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int)
    AVS_GetNrOfDevices = prototype(("AVS_GetNrOfDevices", lib),)
    ret = AVS_GetNrOfDevices()
    return ret

def AVS_GetList(listsize, requiredsize, IDlist):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(AvsIdentityType))
    paramflags = (1, "listsize",), (2, "requiredsize",), (2, "IDlist",),
    AVS_GetList = prototype(("AVS_GetList", lib), paramflags)
    ret = AVS_GetList(listsize)
    # looks like you only pass the '1' parameters here
    # the '2' parameters are returned in 'ret' !!!
    return ret

def AVS_Activate(deviceID):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.POINTER(AvsIdentityType))
    paramflags = (1, "deviceId",),
    AVS_Activate = prototype(("AVS_Activate", lib), paramflags)
    ret = AVS_Activate(deviceID)
    return ret

def AVS_UseHighResAdc(handle, enable):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_bool)
    paramflags = (1, "handle",), (1, "enable",),
    AVS_UseHighResAdc = prototype(("AVS_UseHighResAdc", lib), paramflags)
    ret = AVS_UseHighResAdc(handle, enable)
    return ret

def AVS_GetVersionInfo(handle, FPGAversion, FWversion, DLLversion):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_char * VERSION_LEN, ctypes.c_char * VERSION_LEN, ctypes.c_char * VERSION_LEN)
    paramflags = (1, "handle",), (2, "FPGAversion",), (2, "FWversion",), (2, "DLLversion",),
    AVS_GetVersionInfo =prototype(("AVS_GetVersionInfo", lib), paramflags) 
    ret = AVS_GetVersionInfo(handle)
    return ret

def AVS_PrepareMeasure(handle, measconf):
    lib = ctypes.WinDLL(dll)
    datatype = ctypes.c_byte * 41
    data = datatype()
    temp = datatype()
    temp = struct.pack("HHfIIBBHBBBBBHIIfH", measconf.m_StartPixel,
                                             measconf.m_StopPixel,
                                             measconf.m_IntegrationTime,
                                             measconf.m_IntegrationDelay,
                                             measconf.m_NrAverages,
                                             measconf.m_CorDynDark_m_Enable,
                                             measconf.m_CorDynDark_m_ForgetPercentage,
                                             measconf.m_Smoothing_m_SmoothPix,
                                             measconf.m_Smoothing_m_SmoothModel,
                                             measconf.m_SaturationDetection,
                                             measconf.m_Trigger_m_Mode,
                                             measconf.m_Trigger_m_Source,
                                             measconf.m_Trigger_m_SourceType,
                                             measconf.m_Control_m_StrobeControl,
                                             measconf.m_Control_m_LaserDelay,
                                             measconf.m_Control_m_LaserWidth,
                                             measconf.m_Control_m_LaserWaveLength,
                                             measconf.m_Control_m_StoreToRam )

# copy bytes from temp to data, otherwise you will get a typing error below
# why is this necessary?? they have the same type to start with ??
    x = 0
    while (x < 41): # 0 through 40
        data[x] = temp[x]
        x += 1
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_byte * 41)
    paramflags = (1, "handle",), (1, "measconf",),  
    AVS_PrepareMeasure = prototype(("AVS_PrepareMeasure", lib), paramflags)
    ret = AVS_PrepareMeasure(handle, data)
    return ret

def AVS_Measure(handle, windowhandle, nummeas):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.wintypes.HWND, ctypes.c_uint16)
    paramflags = (1, "handle",), (1, "windowhandle",), (1, "nummeas"),
    AVS_Measure = prototype(("AVS_Measure", lib), paramflags)
    ret = AVS_Measure(handle, windowhandle, nummeas) 
    return ret

class callbackclass(QObject):
    newdata = pyqtSignal()
    def __init__(self):
        QObject.__init__(self, parent)
        self.newdata.connect(PyQt5_demo.MainWindow.handle_newdata)
    def callback(self, handle, error):
        self.newdata.emit() # signal must be from a class !!

# We have not succeeded in getting the callback to execute without problem
# please use AVS_Measure instead using Windows messaging or polling

def AVS_MeasureCallback(handle, adres, nummeas):
    CBTYPE = ctypes.CFUNCTYPE(None, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, CBTYPE, ctypes.c_uint16)
    paramflags = (1, "handle",), (1, "adres",), (1, "nummeas"),
    AVS_MeasureCallback = prototype(("AVS_MeasureCallback", lib), paramflags)
    ret = AVS_MeasureCallback(handle, CBTYPE(callbackclass.callback), nummeas)  # CRASHES python

def AVS_StopMeasure(handle):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int)
    paramflags = (1, "handle",),
    AVS_StopMeasure = prototype(("AVS_StopMeasure", lib), paramflags)
    ret = AVS_StopMeasure(handle)
    return ret

def AVS_PollScan(handle):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_bool, ctypes.c_int)
    paramflags = (1, "handle",),
    AVS_PollScan = prototype(("AVS_PollScan", lib), paramflags)
    ret = AVS_PollScan(handle)
    return ret
    
def AVS_GetScopeData(handle, timelabel, spectrum):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_double * 4096))
    paramflags = (1, "handle",), (2, "timelabel",), (2, "spectrum",),
    AVS_GetScopeData = prototype(("AVS_GetScopeData", lib), paramflags)
    ret = AVS_GetScopeData(handle)
    return ret

def AVS_GetLambda(handle, wavelength):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double * 4096))
    paramflags = (1, "handle",), (2, "wavelength",),
    AVS_GetLambda = prototype(("AVS_GetLambda", lib), paramflags)
    ret = AVS_GetLambda(handle)
    return ret

def AVS_SetDigOut(handle, portId, value):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_uint8, ctypes.c_uint8)
    paramflags = (1, "handle",), (1, "portId",), (1, "value",),
    AVS_SetDigOut = prototype(("AVS_SetDigOut", lib), paramflags)
    ret = AVS_SetDigOut(handle, portId, value)
    return ret

def AVS_GetAnalogIn(handle, AnalogInId, AnalogIn):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_uint8, ctypes.POINTER(ctypes.c_float))
    paramflags = (1, "handle",), (1, "AnalogInId",), (2, "AnalogIn",),
    AVS_GetAnalogIn = prototype(("AVS_GetAnalogIn", lib), paramflags)
    ret = AVS_GetAnalogIn(handle, AnalogInId)
    return ret

def AVS_GetParameter(handle, size, reqsize, deviceconfig):
    lib = ctypes.WinDLL(dll)
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_uint32, ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(DeviceConfigType))
    paramflags = (1, "handle",), (1, "size",), (2, "reqsize",), (2, "deviceconfig",),
    AVS_GetParameter = prototype(("AVS_GetParameter", lib), paramflags)
    ret = AVS_GetParameter(handle, size)
    return ret

def AVS_SetParameter(handle, deviceconfig):
    lib = ctypes.WinDLL(dll)
    datatype = ctypes.c_byte * 63484
    data = datatype()
    temp = datatype()
    temp = struct.pack("HH64B" +
                       "BH5f?8ddd2ff2ff30H" +      # Detector
                       "HBf4096fBI" +              # Irradiance
                       "HBf4096f" +                # Reflectance
                       "4096f" +                   # SpectrumCorrect
                       "?HHfIIBBHBBBBBHIIfHH12B" + # StandAlone 
                       "5f5f5f" +                  # Temperature
                       "?f2f" +                    # TecControl
                       "2f2f10f10f " +             # ProcessControl
                       "IIIBHB" +                  # EthernetSettings
                       "9720B" +                   # Reserved
                       "4096B",                    # OemData
                       deviceconfig.m_Len,
                       deviceconfig.m_ConfigVersion,
                       deviceconfig.m_aUserFriendlyId,
                       deviceconfig.m_Detector_m_SensorType,
                       deviceconfig.m_Detector_m_NrPixels,
                       deviceconfig.m_Detector_m_aFit,
                       deviceconfig.m_Detector_m_NLEnable,
                       deviceconfig.m_Detector_m_aNLCorrect,
                       deviceconfig.m_Detector_m_aLowNLCounts,
                       deviceconfig.m_Detector_m_aHighNLCounts,
                       deviceconfig.m_Detector_m_Gain,
                       deviceconfig.m_Detector_m_Reserved,
                       deviceconfig.m_Detector_m_Offset,
                       deviceconfig.m_Detector_m_ExtOffset,
                       deviceconfig.m_Detector_m_DefectivePixels,
                       deviceconfig.m_Irradiance_m_IntensityCalib_m_Smoothing_m_SmoothPix,
                       deviceconfig.m_Irradiance_m_IntensityCalib_m_Smoothing_m_SmoothModel,
                       deviceconfig.m_Irradiance_m_IntensityCalib_m_CalInttime,
                       deviceconfig.m_Irradiance_m_IntensityCalib_m_aCalibConvers,
                       deviceconfig.m_Irradiance_m_CalibrationType,
                       deviceconfig.m_Irradiance_m_FiberDiameter,
                       deviceconfig.m_Reflectance_m_Smoothing_m_SmoothPix,
                       deviceconfig.m_Reflectance_m_Smoothing_m_SmoothModel,
                       deviceconfig.m_Reflectance_m_CalInttime,
                       deviceconfig.m_Reflectance_m_aCalibConvers,
                       deviceconfig.m_SpectrumCorrect,
                       deviceconfig.m_StandAlone_m_Enable,
                       deviceconfig.m_StandAlone_m_Meas_m_StartPixel,
                       deviceconfig.m_StandAlone_m_Meas_m_StopPixel,
                       deviceconfig.m_StandAlone_m_Meas_m_IntegrationTime,
                       deviceconfig.m_StandAlone_m_Meas_m_IntegrationDelay,
                       deviceconfig.m_StandAlone_m_Meas_m_NrAverages,
                       deviceconfig.m_StandAlone_m_Meas_m_CorDynDark_m_Enable, 
                       deviceconfig.m_StandAlone_m_Meas_m_CorDynDark_m_ForgetPercentage,
                       deviceconfig.m_StandAlone_m_Meas_m_Smoothing_m_SmoothPix,
                       deviceconfig.m_StandAlone_m_Meas_m_Smoothing_m_SmoothModel,
                       deviceconfig.m_StandAlone_m_Meas_m_SaturationDetection,
                       deviceconfig.m_StandAlone_m_Meas_m_Trigger_m_Mode,
                       deviceconfig.m_StandAlone_m_Meas_m_Trigger_m_Source,
                       deviceconfig.m_StandAlone_m_Meas_m_Trigger_m_SourceType,
                       deviceconfig.m_StandAlone_m_Meas_m_Control_m_StrobeControl,
                       deviceconfig.m_StandAlone_m_Meas_m_Control_m_LaserDelay,
                       deviceconfig.m_StandAlone_m_Meas_m_Control_m_LaserWidth,
                       deviceconfig.m_StandAlone_m_Meas_m_Control_m_LaserWaveLength,
                       deviceconfig.m_StandAlone_m_Meas_m_Control_m_StoreToRam,
                       deviceconfig.m_StandAlone_m_Nmsr,
                       deviceconfig.m_StandAlone_m_Reserved,
                       deviceconfig.m_Temperature_1_m_aFit,
                       deviceconfig.m_Temperature_2_m_aFit,
                       deviceconfig.m_Temperature_3_m_aFit,
                       deviceconfig.m_TecControl_m_Enable,
                       deviceconfig.m_TecControl_m_Setpoint,
                       deviceconfig.m_TecControl_m_aFit,
                       deviceconfig.m_ProcessControl_m_AnalogLow,
                       deviceconfig.m_ProcessControl_m_AnalogHigh,
                       deviceconfig.m_ProcessControl_m_DigitalLow,
                       deviceconfig.m_ProcessControl_m_DigitalHigh,
                       deviceconfig.m_EthernetSettings_m_IpAddr,
                       deviceconfig.m_EthernetSettings_m_NetMask,
                       deviceconfig.m_EthernetSettings_m_Gateway,
                       deviceconfig.m_EthernetSettings_m_DhcpEnabled,
                       deviceconfig.m_EthernetSettings_m_TcpPort,
                       deviceconfig.m_EthernetSettings_m_LinkStatus,
                       deviceconfig.m_Reserved,
                       deviceconfig.m_OemData)
    x = 0
    while (x < 63484): # 0 through 63483
        data[x] = temp[x]
        x += 1
    prototype = ctypes.WINFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_byte * 63484)
    paramflags = (1, "handle",), (1, "deviceconfig",),  
    AVS_SetParameter = prototype(("AVS_SetParameter", lib), paramflags)
    ret = AVS_SetParameter(handle, data)
    return ret

if __name__ == '__main__':
        ret = AVS_Init(0)    
        # QMessageBox.information(self,"Info","AVS_Init returned:  {0:d}".format(ret))
        ret = AVS_GetNrOfDevices()
        # QMessageBox.information(self,"Info","AVS_GetNrOfDevices returned:  {0:d}".format(ret))
        req = 0
        mylist = AvsIdentityType * 1
        ret = AVS_GetList(75, req, mylist)
        serienummer = str(ret[1].SerialNumber.decode("utf-8"))
        print("Found Serialnumber: " + serienummer)
        dev_handle = AVS_Activate(ret[1])

        devcon = DeviceConfigType
        reqsize = 0
        ret = AVS_GetParameter(dev_handle, 63484, reqsize, devcon)
        print(ret)
        pixels = ret[1].m_Detector_m_NrPixels
        print(pixels)
        ret = AVS_GetLambda(dev_handle,globals.wavelength)
        x = 0
        while (x < pixels): # 0 through 2047
            globals.wavelength[x] = ret[x]
            x += 1
        #print(globals.wavelength)

        ret = AVS_UseHighResAdc(dev_handle, True)
        measconfig = MeasConfigType
        measconfig.m_StartPixel = 0
        measconfig.m_StopPixel = pixels - 1
        measconfig.m_IntegrationTime = 5
        measconfig.m_IntegrationDelay = 0
        measconfig.m_NrAverages = 100
        measconfig.m_CorDynDark_m_Enable = 0  # nesting of types does NOT work!!
        measconfig.m_CorDynDark_m_ForgetPercentage = 0
        measconfig.m_Smoothing_m_SmoothPix = 0
        measconfig.m_Smoothing_m_SmoothModel = 0
        measconfig.m_SaturationDetection = 0
        measconfig.m_Trigger_m_Mode = 0
        measconfig.m_Trigger_m_Source = 0
        measconfig.m_Trigger_m_SourceType = 0
        measconfig.m_Control_m_StrobeControl = 0
        measconfig.m_Control_m_LaserDelay = 0
        measconfig.m_Control_m_LaserWidth = 0
        measconfig.m_Control_m_LaserWaveLength = 0.0
        measconfig.m_Control_m_StoreToRam = 0
        ret = AVS_PrepareMeasure(dev_handle, measconfig)
        
        # to use Windows messages, supply a window handle to send the messages to
        # ret = AVS_Measure(globals.dev_handle, int(self.winId()), nummeas)
        # single message sent from DLL, confirmed with Spy++

        # when using polling, just pass a 0 for the windows handle
        import time
        import matplotlib.pyplot as plt
        import csv
        import numpy as np
        
        def measure():
            ret = AVS_Measure(dev_handle, 0, 1)
            dataready = False
            while (dataready == False):
                dataready = (AVS_PollScan(dev_handle) == True)
                time.sleep(0.001)
            if dataready == True:
                timestamp = 0
                ret = AVS_GetScopeData(dev_handle, timestamp, globals.spectraldata )
                timestamp = ret[0]
                x = 0
                while (x < pixels): # 0 through 2047
                    globals.spectraldata[x] = ret[1][x]
                    x += 1
            return np.asarray([globals.wavelength[0:pixels], globals.spectraldata[0:pixels]])

        ref = measure()
        # plt.plot(ref[0], ref[1])
        # plt.show()

        nummeas = 30
        data = []
        for i in range(nummeas):
            time.sleep(10)
            d = measure()
            d[1] = d[1]/ref[1]
            data.append(d)
            # plt.plot(data[i][0], data[i][1])
            # plt.show()

        fname = 'C:/Users/MatterLab/Dropbox/PythonLab/pylab/instruments/Avaspec/test_DHmini_direct_expsure5_average100_interval10s.csv'
        out = np.asarray([data[i][1] for i in range(nummeas)])
        print(ref.shape, out.shape)
        out = np.transpose(np.concatenate((ref, out)))
        head = ['wavelength', 'ref']
        head.extend(['%s s' %(i*10) for i in range(nummeas)])
         
        with open(fname, 'w', newline ='') as f:
            writer = csv.writer(f)
            writer.writerow(head)
            writer.writerows(out)   