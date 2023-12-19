import pylab.instruments as ins
from pylab.instruments import pv_const
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
from statistics import mean


def set_filter():   
    ps._startSignal()
    time.sleep(2)
    while True:
        pos = fw.get_position()
        s_count, c_count = th260.get_sync_rate(), th260.get_count_rate(0)
        print(pos, s_count, c_count)

        if th261.get_count_rate(0) <  1e-2 * th260.get_sync_rate():
           if pos == 2:
               print('Warning : laser power is too weak')
               return 1
           else:
                fw.set_position(pos-1)
        else: 
            if pos == 12:
                print('Warning : laser power is too strong')
                return 1
            else:
                fw.set_position(pos+1)
                break
    ps._stopSignal()


def exp_fit(x, hist, rise = 50):
    def exp(x, amp, tau, x0):
        return (amp * np.exp(-1*(x-x0)/tau))**fit_order
    c_max, max_index = np.max(hist), np.argmax(hist)
    bg = np.average(hist[0:max_index-rise])
    a_hist = hist[max_index:]-bg
    a_hist[a_hist<0]=0.00001   
    a_x, a_hist = x[max_index:]-x[max_index], (a_hist)**fit_order
    init_vals = [c_max, 10, 0]
    fit_vals, covar = curve_fit(exp, a_x, a_hist, p0=init_vals, bounds=(0,[np.inf, np.inf, 4]))
    print(fit_vals, bg)
    return fit_vals[0], fit_vals[1], fit_vals[2], bg, max_index

def exp_fit_log(x, hist, rise = 50):
    def exp_log(x, amp, tau, x0):
        return np.log(amp * np.exp(-1*(x-x0)/tau))
    c_max, max_index = np.max(hist), np.argmax(hist)
    bg = np.average(hist[0:max_index-rise])
    a_hist = hist[max_index:]-bg
    a_hist[a_hist<0]=0.00001
    a_x, a_hist = x[max_index:]-x[max_index], np.log(a_hist)
    plt.plot(a_x,a_hist)
    plt.show()
    init_vals = [c_max, 10, 0]
    fit_vals, covar = curve_fit(exp_log, a_x, a_hist, p0=init_vals, bounds=(0,[np.inf, np.inf, 2]))
    print(fit_vals, bg)
    return fit_vals[0], fit_vals[1], fit_vals[2], bg, max_index

def measure(overflow, timeout):
    th260.set_overflow(overflow)
    ps._startSignal()
    time.sleep(1)
    print(th260.get_sync_rate())
    hist = np.asarray(th260.measure(timeout))
    ps._stopSignal()
    x = np.asarray([1e-3 * i * res * 2**binning  for i in range(len(hist))])
    return x, hist

fit_order = 1

fw = ins.ThorlabsFW212C('ASRL2::INSTR', pos_count=12)
fw.set_position(12)
##
# initialize picoscope
ps = ins.PS5242D("8BIT")
ps._setBuiltInSignal(1, 1e6, 1500, 3000)

# initialize timeharp
th260 = ins.TH260(0) 
print(th260.get_lib_version())
print(th260.get_hardware_info())
print(th260.get_numchannels())

th260.set_sync_divider(1)
th260.set_cfd(channel='SYNC', cfd_level=-70, zero_cross=-10)
th260.set_cfd(channel=0,      cfd_level=-70, zero_cross=-10)

th260.set_channel_offset(channel='SYNC', offset=0)
th260.set_channel_offset(channel=0,      offset=0)

binning = 0
th260.set_binning(binning)
th260.set_offset(0)

res = th260.get_resolution()

#set_filter
if set_filter() == 1:
    exit()

#set_frequency 
x, hist = measure(100, 60)
amp, tau, x0, bg, max_index = exp_fit(x, hist)
x_fit = x[max_index:]
hist_fit = [amp*np.exp(-1*(xx-x[max_index]-x0)/tau)+bg for xx in x_fit]

frequency = int(min(1/(tau*9.212*3*1e-9), 2e7))
print(frequency)
ps._setBuiltInSignal(1, frequency, 1500, 3000)

#Measure Histgram 
x, hist = measure(10000, 60)
amp, tau, x0, bg, max_index = exp_fit(x, hist)
x_fit = x[max_index:]
hist_fit = [amp*np.exp(-1*(xx-x[max_index]-x0)/tau)+bg for xx in x_fit]


#plot results
plt.plot(x, hist)
plt.plot(x_fit, hist_fit)
plt.xlim(10,0.9*1e9/frequency)
plt.ylim(1e0,2e4)
plt.yscale('log')
plt.xlabel('time/ns')
plt.ylabel('counts')
plt.text(15,5e3, 'tau = {} ns'.format(tau))
plt.show()

th260.close_all()
ps._close()





