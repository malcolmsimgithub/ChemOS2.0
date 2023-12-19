import pylab.instruments as ins
from pylab.instruments import pv_const
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time
from statistics import mean


def set_filter():   
    while True:
        pos = fw.get_position()
        s_count,c_count = th260.get_sync_rate(),th260.get_count_rate(0)
        print(pos, s_count, c_count)

        if th260.get_count_rate(0) <  1e-2 * th260.get_sync_rate():
           if pos == 1:
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

def exp_fit(x, hist, rise = 50):
    def exp(x, amp, tau):
        return amp * np.exp(-1*x/tau)
    c_max, max_index = np.max(hist), np.argmax(hist)
    bg = np.average(hist[0:max_index-rise])
    a_x, a_hist = x[max_index:]-x[max_index], hist[max_index:]-bg
    init_vals = [c_max, 10]
    fit_vals, covar = curve_fit(exp, a_x, a_hist, p0=init_vals)
    print(fit_vals, bg)
    return fit_vals[0], fit_vals[1], bg, max_index

def exp_fit_deconv(org_x, hist, rise = 50):
    def conv_exp(x, amp, tau):
        flt = np.loadtxt('./IRF_th260P.txt')[1]
        return np.convolve(amp * np.exp(-1*x/tau), flt, mode = 'same')
    c_max, max_index = np.max(hist), np.argmax(hist)
    bg = np.average(hist[0:max_index-rise])
    a_x, a_hist = x[max_index:]-x[max_index], hist[max_index:]-bg
    init_vals = [c_max, 1]
    fit_vals, covar = curve_fit(conv_exp, a_x, a_hist, p0=init_vals)

    print(fit_vals, bg)
    return fit_vals[0], fit_vals[1], bg, max_index

fw = ins.ThorlabsFW212C('ASRL2::INSTR', pos_count=12)
fw.set_position(12)
##
# initialize timeharp

th260 = ins.TH260(0)
print(th260.get_lib_version())
print(th260.get_hardware_info())
print(th260.get_numchannels())

th260.set_sync_divider(1)
th260.set_cfd(channel='SYNC', cfd_level=-50, zero_cross=-10)
th260.set_cfd(channel=0,      cfd_level=-50, zero_cross=-10)

th260.set_channel_offset(channel='SYNC', offset=0)
th260.set_channel_offset(channel=0,      offset=0)

binning = 0
th260.set_binning(binning)
th260.set_offset(0)
th260.set_overflow(10000)

res = th260.get_resolution()

#set filter
if set_filter() == 1:
    exit()

#Measure histgram
hist = np.asarray(th260.measure(30))
x = np.asarray([1e-3 * i * res * 2**binning  for i in range(len(hist))])

#exp_fit
amp, tau, bg, max_index = exp_fit_deconv(x, hist)
x_fit = x[max_index:]
flt = np.loadtxt('./IRF_th260P.txt')[1]
hist_fit = [amp * np.exp(-1*(xx-x[max_index])/tau)+bg for xx in x_fit]
hist_fit_conv = np.convolve(hist_fit, flt, mode='same')

#plot results
plt.plot(x, hist)
plt.plot(x_fit, hist_fit)
plt.plot(x_fit, hist_fit_conv)
plt.xlim(10,100)
plt.ylim(1e0,1e5)
plt.yscale('log')
plt.xlabel('time/ns')
plt.ylabel('counts')
plt.text(15,5e4, 'tau = {} ns'.format(tau))
plt.show()


th260.close_all()






