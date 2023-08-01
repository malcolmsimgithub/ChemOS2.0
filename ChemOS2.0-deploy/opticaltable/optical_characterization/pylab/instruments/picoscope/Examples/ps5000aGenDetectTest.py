#
# Copyright (C) 2018 Pico Technology Ltd. See LICENSE file for terms.
#
# PicoScope 5000 (A API) Signal Generator Example
# This example demonstrates how to use the PicoScope 5000 Series (ps5000a) driver API functions to set up the signal generator to do the following:
# 
# 1. Output a sine wave 
# 2. Output a square wave 
# 3. Output a sweep of a square wave signal

import ctypes
from picosdk.ps5000a import ps5000a as ps
import time
import numpy as np
import matplotlib.pyplot as plt
from picosdk.functions import adc2mV, assert_pico_ok, mV2adc

status = {}
chandle = ctypes.c_int16()

resolution =ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_8BIT"]

# Open the device
status["openUnit"] = ps.ps5000aOpenUnit(ctypes.byref(chandle), None, resolution)


try:
    assert_pico_ok(status["openUnit"])
except: # PicoNotOkError:

    powerStatus = status["openUnit"]

    if powerStatus == 286:
        status["changePowerSource"] = ps.ps5000aChangePowerSource(chandle, powerStatus)
    elif powerStatus == 282:
        status["changePowerSource"] = ps.ps5000aChangePowerSource(chandle, powerStatus)
    else:
        raise

    assert_pico_ok(status["changePowerSource"])


# Set up channel A
# handle = chandle
channel = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
# enabled = 1
coupling_type = ps.PS5000A_COUPLING["PS5000A_DC"]
chARange = ps.PS5000A_RANGE["PS5000A_5V"]
# analogue offset = 0 V
status["setChA"] = ps.ps5000aSetChannel(chandle, channel, 1, coupling_type, chARange, 0)
assert_pico_ok(status["setChA"])

maxADC = ctypes.c_int16()
status["maximumValue"] = ps.ps5000aMaximumValue(chandle, ctypes.byref(maxADC))
assert_pico_ok(status["maximumValue"])

# Set up single trigger
# handle = chandle
# enabled = 1
source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
threshold = int(mV2adc(500,chARange, maxADC))
print(threshold)
# direction = PS5000A_RISING = 2
# delay = 0 s
# auto Trigger = 1000 ms
status["trigger"] = ps.ps5000aSetSimpleTrigger(chandle, 1, source, threshold, 2, 0, 1000)
assert_pico_ok(status["trigger"])

# Set number of pre and post trigger samples to be collected
preTriggerSamples = 2500
postTriggerSamples = 2500
maxSamples = preTriggerSamples + postTriggerSamples

# Get timebase information
# handle = chandle
#----------------------------------------------------------------------
#sampling rate
#resolution : 8bit ->  timebase 1:1ns, 2:2ns, 3:4ns, 4:8ns, ....
#            12bit ->  timebase 1:2ns, 2:4ns, 3:8ns, 4:16ns, ...
#            14bit,15bit ->  timebase 3: 8ns, 4:16ns
#            16bit t> timebase 4:16ns, 8: 32ns.....
#----------------------------------------------------------------------

timebase = 1  
# noSamples = maxSamples
# pointer to timeIntervalNanoseconds = ctypes.byref(timeIntervalns)
# pointer to maxSamples = ctypes.byref(returnedMaxSamples)
# segment index = 0
timeIntervalns = ctypes.c_float()
returnedMaxSamples = ctypes.c_int32()
status["getTimebase2"] = ps.ps5000aGetTimebase2(chandle, timebase, maxSamples, ctypes.byref(timeIntervalns), ctypes.byref(returnedMaxSamples), 0)
assert_pico_ok(status["getTimebase2"])

#Generator setting

# 0:sine, 1:square, 2:triangle, 3:Ramp_Up, 4:Ramp_down, 5:sinc, 6:gausian, 7:h_sine, 8:DC, 9:whitenoise
wavetype = ctypes.c_int32(1) 
sweepType = ctypes.c_int32(0) #0:UP, 1:Down, 2:UPDOWN, 3:DOWNUP, 4:Max_Sweep_type
triggertype = ctypes.c_int32(0) # 0:Rising, 1:Falling, 2:Gate_high, 3:Gate_low
triggerSource = ctypes.c_int32(0) # 0:None, 1:Scope_trig, 2:Aux_in, 3:Ext_in, 4:Software_trig
offset = 1000000 #uv
pkTopk = 2000000 #uV
startFreq = 5e6 #Hz
stopFreq = 5e6 #Hz
increment = 0  #increment of frequency
dwelltime = 1 #s, time to stay at each frequency
operation = 0 #0:Normal, 1:whitenoise, 
shots = 0 #0:sweep as specified by sweeps, 1: maximum number
sweeps = 0 # 0:produce cycles specified by shot, 1: maximum number
extInThreshold = 0 #Ext_trigger level

status["setSigGenBuiltInV2"] = ps.ps5000aSetSigGenBuiltInV2(chandle, offset, pkTopk, wavetype, startFreq, stopFreq, increment, dwelltime, sweepType, operation, shots, sweeps, triggertype, triggerSource, extInThreshold)
assert_pico_ok(status["setSigGenBuiltInV2"])

#run block detection
status["runBlock"] = ps.ps5000aRunBlock(chandle, preTriggerSamples, postTriggerSamples, timebase, None, 0, None, None)
assert_pico_ok(status["runBlock"])

# Check for data collection to finish using ps5000aIsReady
ready = ctypes.c_int16(0)
check = ctypes.c_int16(0)
while ready.value == check.value:
    status["isReady"] = ps.ps5000aIsReady(chandle, ctypes.byref(ready))

# Create buffers ready for assigning pointers for data collection
bufferAMax = (ctypes.c_int16 * maxSamples)()
bufferAMin = (ctypes.c_int16 * maxSamples)() # used for downsampling which isn't in the scope of this example

# Set data buffer location for data collection from channel A
# handle = chandle
source = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
# pointer to buffer max = ctypes.byref(bufferAMax)
# pointer to buffer min = ctypes.byref(bufferAMin)
# buffer length = maxSamples
# segment index = 0
# ratio mode = PS5000A_RATIO_MODE_NONE = 0
status["setDataBuffersA"] = ps.ps5000aSetDataBuffers(chandle, source, ctypes.byref(bufferAMax), ctypes.byref(bufferAMin), maxSamples, 0, 0)
assert_pico_ok(status["setDataBuffersA"])

# create overflow loaction
overflow = ctypes.c_int16()
# create converted type maxSamples
cmaxSamples = ctypes.c_int32(maxSamples)

# Retried data from scope to buffers assigned above
# handle = chandle
# start index = 0
# pointer to number of samples = ctypes.byref(cmaxSamples)
# downsample ratio = 0
# downsample ratio mode = PS5000A_RATIO_MODE_NONE
# pointer to overflow = ctypes.byref(overflow))
status["getValues"] = ps.ps5000aGetValues(chandle, 0, ctypes.byref(cmaxSamples), 0, 0, 0, ctypes.byref(overflow))
assert_pico_ok(status["getValues"])


# convert ADC counts data to mV
adc2mVChAMax =  adc2mV(bufferAMax, chARange, maxADC)

# Create time data
time = np.linspace(0, (cmaxSamples.value) * timeIntervalns.value, cmaxSamples.value)

# plot data from channel A and B
plt.plot(time, adc2mVChAMax[:])
plt.xlabel('Time (ns)')
plt.ylabel('Voltage (mV)')
plt.show()

# Stop the scope
# handle = chandle
status["stop"] = ps.ps5000aStop(chandle)
assert_pico_ok(status["stop"])


#time.sleep(10)


# Closes the unit
# Handle = chandle
status["close"] = ps.ps5000aCloseUnit(chandle)
#status["stop"] = ps.ps5000aStop(chandle)
#assert_pico_ok(status["stop"])

# Displays the status returns
print(status)
