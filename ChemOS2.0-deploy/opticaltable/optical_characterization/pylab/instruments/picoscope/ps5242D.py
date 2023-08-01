import time, os
import numpy as np
from .picosdk.ps5000a2 import ps5000a
from .picosdk.functions import adc2mV, assert_pico_ok, mV2adc
import matplotlib.pyplot as plt
import math

import ctypes

class PS5242D(ps5000a):
    """ 
    This is a Python module defining the functions from the ps5000aApi.h C header
    file for PicoScope 5000 Series oscilloscopes using the ps5000a driver API
    functions. 
    """
    def __init__(self, resolution = "8BIT"):
        
        super(PS5242D, self).__init__()
        self.status = {}
        self.chandle = ctypes.c_int16()
        self.resolution = self.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_" + resolution]
        self.chARange = self.PS5000A_RANGE["PS5000A_20V"]
        self.chBRange = self.PS5000A_RANGE["PS5000A_20V"]
        self.maxADC = ctypes.c_int16()
        self.status["openunit"] = self.ps5000aOpenUnit(ctypes.byref(self.chandle), None, self.resolution)
        try:
            assert_pico_ok(self.status["openunit"])
        except: # PicoNotOkError:
            powerStatus = self.status["openunit"]

            if powerStatus == 286:
                self.status["changePowerSource"] = self.ps5000aChangePowerSource(self.chandle, powerStatus)
            elif powerStatus == 282:
                self.status["changePowerSource"] = self.ps5000aChangePowerSource(self.chandle, powerStatus)
            else:
                raise
            assert_pico_ok(self.status["changePowerSource"])

        self.status["maximumValue"] = self.ps5000aMaximumValue(self.chandle, ctypes.byref(self.maxADC))
        assert_pico_ok(self.status["maximumValue"])


    def setResolution(self, resolution):
        """  
        Set resolution of the scope

        Parameters
        --------------------
        resolution : str
            '8BIT', '12BIT', '14BIT', '15BIT',  '16BIT'
        """
        self.resolution = self.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_" + resolution]
        self.ps5000aSetDeviceResolution(self.chandle, self.resolution)
#        self.res = ctypes.c_int16()
#        self.ps5000aGetDeviceResolution(self.chandle, ctypes.byref(self.res))
#        print(self.res)
#        self.status["openunit"] = self.ps5000aOpenUnit(ctypes.byref(self.chandle), None, self.resolution)


    def _calcTimebase(self, interval, resolution):
        """  
        Calculate the timebase of the scope

        Parameters
        --------------------
        interval : float
            sampling interval in microseconds
        resolution : str
            '8BIT', '12BIT', '14BIT', '15BIT',  '16BIT'

        Reterns
        --------------------
        timebase : int
        interval : float
        """
        #interval in us
        if resolution == "8BIT":
            if interval < 0.008:
                timebase = int(math.log(interval*1000, 2))
                interval = 0.001* 2**timebase
            else:
                timebase = int(interval*1000/8+2)
                interval = 8*(timebase-2)/1000
        if resolution == "12BIT":
            if interval < 0.016:
                timebase = int(math.log(interval*1000, 2))
                interval = 0.001 * 2**timebase
            else:
                timebase = int(interval*1000/16+3)
                interval = 16*(timebase-3)/1000 
        if resolution == "14BIT":
            if interval < 0.016:
                timebase = int(math.log(interval*1000, 2))
                interval = 0.001 * 2**timebase
            else:
                timebase = int(interval*1000/8+2)
                interval = 16*(timebase-2)/1000
        return timebase, interval


    def setChannel(self, channel = "A", coupling = "DC", Crange = "20V"):
        """  
        Setup channels of the scope.

        Parameters
        --------------------
        channel : str
            channel to set ("A", "B", ..)
        coupling : str
            'DC' or 'AC' 
        Crange : str
            '10MV', '20MV', '50MV', '100MV', '200MV', '500MV', '1V', '2V', '5V', '10V', '20V' 
        """
        chl = self.PS5000A_CHANNEL["PS5000A_CHANNEL_" + channel]
        coupling_type = self.PS5000A_COUPLING["PS5000A_" + coupling]

        if channel == "A":
            self.chARange = self.PS5000A_RANGE["PS5000A_" + Crange]
            # analogue offset = 0 V
            self.status["setChA"] = self.ps5000aSetChannel(self.chandle, chl, 1, coupling_type, self.chARange, 0)
            assert_pico_ok(self.status["setChA"])
        elif channel == "B":
            self.chBRange = self.PS5000A_RANGE["PS5000A_" + Crange]
            self.status["setChB"] = self.ps5000aSetChannel(self.chandle, chl, 1, coupling_type, self.chBRange, 0)
            assert_pico_ok(self.status["setChB"])
        print('Picoscope : Channel %s was set to %s%s range' %(channel, coupling, Crange))
            


    def setTrigger(self, channel = "A", threshold = 50, direction = 2, delay = 0, autotrig = 1000):
        """  
        enable the trigger of the scope.

        Parameters
        --------------------
        channel : str
            channel to set ("A", "B", ..)
        threshold : float
            trigger threshold in mV
        direction : int
            0: above, 1: below, 2: rising, 3: falling, 4: rising or falling,...
        delay : int
            number of the sample scope would wait before sampling
        autotrig : int
            the number of milliseconds after which the device starts capturing if no trigger occurs.
            If this is set to zero, the scope device waits indefinitely for a trigger
        """
        source = self.PS5000A_CHANNEL["PS5000A_CHANNEL_" + channel]
        if channel == "A":
            crange = self.chARange
        elif channel == "B":
            crange = self.chBRange
        threshold = int(mV2adc(threshold, crange, self.maxADC))

        self.status["trigger"] = self.ps5000aSetSimpleTrigger(self.chandle, 1, 
                                                source, threshold, direction, delay, autotrig)
        assert_pico_ok(self.status["trigger"])



    def disableTrigger(self, channel = 'A'):
        """  
        disable the trigger of the scope.

        Parameters
        --------------------
        channel : str
            channel to set ("A", "B", ..)
        """
        source = self.PS5000A_CHANNEL["PS5000A_CHANNEL_" + channel]
        self.status["trigger"] = self.ps5000aSetSimpleTrigger(self.chandle, 0, 
                                                    source, 0, 2, 0, 1000)
        assert_pico_ok(self.status["trigger"])



    def _runBlock(self, timebase = 1, preTriggerSamples = 2000, postTriggerSamples = 2000):
        """  
        Run scope and collect block of data

        Parameters
        --------------------
        timebase : int
            sampling rate :  
                8bit ->  timebase 1:1ns, 2:2ns, 3:4ns, 4:8ns, ....
                12bit ->  timebase 1:2ns, 2:4ns, 3:8ns, 4:16ns, ...
                14bit,15bit ->  timebase 3: 8ns, 4:16ns
                16bit t> timebase 4:16ns, 8: 32ns.....
        preTriggerSamples : int
            number of sample collected before the trigger events
        postTriggerSamples : int
            number of sample collected after the trigger events

        Reterns
        --------------------
        data : numpy array
            [time, data(channelA), data(channelB)]
        """
        maxSamples = preTriggerSamples + postTriggerSamples

        timeIntervalns = ctypes.c_float()   
        returnedMaxSamples = ctypes.c_int32()
        self.status["getTimebase2"] = self.ps5000aGetTimebase2(self.chandle, timebase, 
                       maxSamples, ctypes.byref(timeIntervalns), ctypes.byref(returnedMaxSamples), 0)
        assert_pico_ok(self.status["getTimebase2"])

        self.status["runBlock"] = self.ps5000aRunBlock(self.chandle, preTriggerSamples, 
                                                    postTriggerSamples, timebase, None, 0, None, None)
        assert_pico_ok(self.status["runBlock"])

        # Check for data collection to finish using ps5000aIsReady
        ready = ctypes.c_int16(0)
        check = ctypes.c_int16(0)
        while ready.value == check.value:
            self.status["isReady"] = self.ps5000aIsReady(self.chandle, ctypes.byref(ready))

        bufferAMax = (ctypes.c_int16 * maxSamples)()
        bufferAMin = (ctypes.c_int16 * maxSamples)()
        source = self.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        self.status["setDataBuffersA"] = self.ps5000aSetDataBuffers(self.chandle, source, 
                                 ctypes.byref(bufferAMax), ctypes.byref(bufferAMin), maxSamples, 0, 0)
        assert_pico_ok(self.status["setDataBuffersA"])

        bufferBMax = (ctypes.c_int16 * maxSamples)()       
        bufferBMin = (ctypes.c_int16 * maxSamples)()
        source = self.PS5000A_CHANNEL["PS5000A_CHANNEL_B"]
        self.status["setDataBuffersB"] = self.ps5000aSetDataBuffers(self.chandle, source, 
                                 ctypes.byref(bufferBMax), ctypes.byref(bufferBMin), maxSamples, 0, 0)
        assert_pico_ok(self.status["setDataBuffersB"])

        overflow = ctypes.c_int16()
        cmaxSamples = ctypes.c_int32(maxSamples)

        self.status["getValues"] = self.ps5000aGetValues(self.chandle, 0, 
                                            ctypes.byref(cmaxSamples), 0, 0, 0, ctypes.byref(overflow))
        assert_pico_ok(self.status["getValues"])

        adc2mVChAMax =  adc2mV(bufferAMax, self.chARange, self.maxADC)
        adc2mVChBMax =  adc2mV(bufferBMax, self.chBRange, self.maxADC)
        
        time = np.linspace(0, (cmaxSamples.value) * timeIntervalns.value, cmaxSamples.value)

        self.status["stop"] = self.ps5000aStop(self.chandle)
        assert_pico_ok(self.status["stop"])

        return np.array([time,adc2mVChAMax,adc2mVChBMax])


    def measure(self, Arange, Brange, resolution, interval, duration, \
                        t_position = 0.5, trigger = None, t_thred = 10, t_dir = 2):
        """  
        Run scope and collect data

        Parameters
        --------------------
        Arange : str
            range of channel A (10MV,20MV,50MV,100MV,200MV,500MV,1V,2V,5V,10V,20V)
        Brange : str
            range of channel B (10MV,20MV,50MV,100MV,200MV,500MV,1V,2V,5V,10V,20V)
        resolution : str
            '8BIT', '12BIT', '14BIT', '15BIT',  '16BIT'
        interval : float
            sampling interval in microseconds
        duration : float
            duration of the data collection in microseconds
        trigger_position : float
            relative trigger position in measurement period (0-1)
        trigger : str or None
            source of the trigger ('A', 'B' or None)
        t_thred : float
            trigger threshold in millivolt 
        t_dir : int
            trigger direction(0: above, 1: below, 2: rising, 3: falling, 4: rising or falling)

        Reterns
        --------------------
        data : numpy array
            [time, data(channelA), data(channelB)]  
            data is in mV
        """
        tbase, interval = self._calcTimebase(interval, resolution)
        self.setResolution(resolution)
        n_sample = int(duration/interval)

        self.setChannel(channel = "A", coupling = "DC", Crange = Arange)
        self.setChannel(channel = "B", coupling = "DC", Crange = Brange)

        if trigger == None:
            self.disableTrigger()
        elif trigger == 'A':
            self.setTrigger(channel = "A", threshold = t_thred, direction = 2, delay = 0) 
        elif trigger == 'B':
            self.setTrigger(channel = "B", threshold = t_thred, direction = 2, delay = 0)
 
        data = self._runBlock(timebase = tbase, preTriggerSamples = int(n_sample*t_position), \
                        postTriggerSamples = int(n_sample*(1-t_position)))

        return data

    def setBuiltInSignal(self, wavetype, frequency, offset, amplitude, sweeptype = 0, triggertype = 2, triggerSource = 4):
        """  
        Setup build-in output signals.

        Parameters
        --------------------
        wavetype : int
            0: sine, 1: square, 2: triangle, 3: ramp_up, 4: ramp_down, 5: sinc, 6: gaussian, 7: half sine, 8: DC, 9: white noise
        frequency : float
            signal frequency in Hz
        offset : float
            siganl offset in millivolts
        amplitude : float
            siganl amplitude in millivolts
        sweeptype : int
            direction of the frequency sweep (0: up, 1: down, 2: updown, 3: downup)
        triggertype : int
            0: rising, 1: falling, 2: gate high, 3: gate low
            signal generator can be started by a rising or falling edge on the trigger signal 
            or can be gated to run whilethe trigger signal is high or low.
        triggerSource : int      
            0: None, 1: scope trigger, 2: aux_in, 3: ext_in, 4: software_trigger
        """
        wavetype = ctypes.c_int32(wavetype)
        sweeptype = ctypes.c_int32(sweeptype)
        triggertype = ctypes.c_int32(triggertype)
        triggerSource = ctypes.c_int32(triggerSource)


        self.status["setSigGenBuiltInV2"] = self.ps5000aSetSigGenBuiltInV2(self.chandle, offset*1000, 
                 amplitude*1000, wavetype, frequency, frequency, 0, 1, sweeptype, 0, 0, 0, triggertype, triggerSource, 0)
        assert_pico_ok(self.status["setSigGenBuiltInV2"])


    def setArbitrarySignal(self, wave_form, frequency, offset, amplitude, 
                                 indexmode = 0, sweeptype = 0, triggertype = 2, triggerSource = 4):
        """  
        Setup arbitrary output signals.

        Parameters
        --------------------
        wave_form : list
            list to discribe the waveform
        frequency : float
            signal frequency in Hz
        offset : float
            signal offset in millivolts
        amplitude : float
            signal peak-to-peak amplitude in millivolts
        indexmode : int
            0: single, 1: dual (generator output the reverse signal after each waveform)
        sweeptype : int
            direction of the frequency sweep (0: up, 1: down, 2: updown, 3: downup)
        triggertype : int
            0: rising, 1: falling, 2: gate high, 3: gate low
            signal generator can be started by a rising or falling edge on the trigger signal 
            or can be gated to run whilethe trigger signal is high or low.
        triggerSource : int      
            0: None, 1: scope trigger, 2: aux_in, 3: ext_in, 4: software_trigger
        """
        MinAmp, MaxAmp, MinSize, MaxSize = self._getAWGMinMax()
        AW_size = len(wave_form)

        if AW_size < MinSize or MaxSize < AW_size:
            raise ValueError('length of the waveform should be {} to {}'.format(MinSize, MaxSize))
            
        frequency = ctypes.c_double(frequency)
        sweeptype = ctypes.c_int32(sweeptype)
        triggertype = ctypes.c_int32(triggertype)
        triggerSource = ctypes.c_int32(triggerSource)
        AW_buffer = (ctypes.c_int16*AW_size)()
        
        print((MaxAmp-MinAmp),(max(wave_form)-min(wave_form)))

        for i in range(AW_size):
            AW_buffer[i] = int((wave_form[i]-min(wave_form))*((MaxAmp-MinAmp)
                                          /(max(wave_form)-min(wave_form)))-((MaxAmp-MinAmp)/2))
        
        phase = self._getAWGPhase(frequency, indexmode, AW_size)
        AW_size = ctypes.c_int32(AW_size)
        self.status["setSigGenArbitrary"] = self.ps5000aSetSigGenArbitrary(self.chandle, offset*1000, amplitude*1000, 
             phase, phase, 0, 1, ctypes.byref(AW_buffer), AW_size, sweeptype, 0, 0, 0, 0, triggertype, triggerSource, 0)


    def _getAWGPhase(self,frequency, indexmode, bufferlength):
        
        phase = ctypes.c_int(32)
        self.ps5000aSigGenFrequencyToPhase(self.chandle, frequency, indexmode, bufferlength, ctypes.byref(phase))

        return phase.value


    def _getAWGMinMax(self):
        minAV, maxAV = ctypes.c_int16(), ctypes.c_int16()
        minAS, maxAS = ctypes.c_int32(), ctypes.c_int32()
        self.ps5000aSigGenArbitraryMinMaxValues(self.chandle, ctypes.byref(minAV), 
                                                    ctypes.byref(maxAV),ctypes.byref(minAS), ctypes.byref(maxAS))

        print(minAV.value, maxAV.value, minAS.value, maxAS.value)

        return minAV.value, maxAV.value, minAS.value, maxAS.value
    

    def makeStepFunction(self, datalength, duty, delay):
        """  
        create step function for arbitrary output.

        Parameters
        --------------------
        datalength: int
            length of the data for the waveform
        duty : float
            duty ratio (0-1)
        delay : float 
            delay of the signal (0-1)
    
        Reterns
        --------------------
        wave_form : list
        """
        if duty < 0 or 1 < duty:
            raise ValueError('values of the duty should be {} to {}'.format(0, 1))
        elif delay < 0 or 1 < delay:
            raise ValueError('values of the delay should be {} to {}'.format(0, 1))

        wave_form = [0 for i in range(datalength)]
        for i in range(int(datalength*duty)):
            wave_form[i+int(datalength*delay)] = 1

        return wave_form


    def startSignal(self):
        """  
        Start signal.
        """
        self.ps5000aSigGenSoftwareControl(self.chandle, 1)


    def stopSignal(self):
        """  
        Stop signal.
        """
        self.ps5000aSigGenSoftwareControl(self.chandle, 0)  
 

    def close(self):
        """  
        Close scope.
        """
        self.status["close"]=self.ps5000aCloseUnit(self.chandle)
        assert_pico_ok(self.status["close"])



if __name__ == "__main__":

    ps5242d = PS5242D("8BIT")
    
#    ps5242d.SetChannel("A","DC", "200MV")
#    ps5242d.SetTrigger("A", -50, 2, 0)
#    data = ps5242d.RunBlock(1, 2000, 2000)
#
#    plt.plot(data[0],data[1])
#    plt.show()

#    ps5242d.OutputBuiltInSignal(1, 5e6,1500,3000)
#    time.sleep(10)

    ps5242d.close()




