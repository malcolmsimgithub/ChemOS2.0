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
from picosdk.functions import assert_pico_ok


status = {}
chandle = ctypes.c_int16()

# Open the device
status["openUnit"] = ps.ps5000aOpenUnit(ctypes.byref(chandle), None, 1)


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

# 0:sine, 1:square, 2:triangle, 3:Ramp_Up, 4:Ramp_down, 5:sinc, 6:gausian, 7:h_sine, 8:DC, 9:whitenoise
wavetype = ctypes.c_int32(1) 
sweepType = ctypes.c_int32(0) #0:UP, 1:Down, 2:UPDOWN, 3:DOWNUP, 4:Max_Sweep_type
triggertype = ctypes.c_int32(0) # 0:Rising, 1:Falling, 2:Gate_high, 3:Gate_low
triggerSource = ctypes.c_int32(0) # 0:None, 1:Scope_trig, 2:Aux_in, 3:Ext_in, 4:Software_trig
offset = 1000000 #uv
pkTopk = 2000000 #uV
startFreq = 1e6 #Hz
stopFreq = 1e6 #Hz
increment = 0  #increment of frequency
dwelltime = 1 #s, time to stay at each frequency
operation = 0 #0:Normal, 1:whitenoise, 
shots = 0 #0:sweep as specified by sweeps, 1: maximum number
sweeps = 0 # 0:produce cycles specified by shot, 1: maximum number
extInThreshold = 0 #Ext_trigger level

status["setSigGenBuiltInV2"] = ps.ps5000aSetSigGenBuiltInV2(chandle, offset, pkTopk, wavetype, startFreq, stopFreq, increment, dwelltime, sweepType, operation, shots, sweeps, triggertype, triggerSource, extInThreshold)
assert_pico_ok(status["setSigGenBuiltInV2"])


# Pauses the script to show signal
time.sleep(10)


# Closes the unit
# Handle = chandle
status["close"] = ps.ps5000aCloseUnit(chandle)
#status["stop"] = ps.ps5000aStop(chandle)
#assert_pico_ok(status["stop"])

# Displays the status returns
print(status)
