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
        s_count = th260.get_sync_rate()
        c_count = th260.get_count_rate(0)
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

def extract_data(data, thread, mode = 'both'):
    c_max = np.max(data)
    max_index = np.argmax(data)
    for i in range (max_index):
        if data[max_index - i] < c_max*thread:
            s_index = max_index -i
            break
    for i in range (len(data)-max_index):
        if data[max_index + i] < c_max*thread:
            e_index = max_index + i
            break
    if mode =='both':
        return data[s_index:e_index+1], [s_index, e_index]
    elif mode == 'forward':
        return data[:e_index+1], e_index
    elif mode == 'backward':
        return data[s_index:], s_index

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
hist = np.asarray(th260.measure(10))

#calculate FWHM
e_data, index =  extract_data(hist, 0.5)
FWHM = (index[1]-index[0])* res * 2**binning

#export IRF
threshold = 0.001
IRF, index= extract_data(hist, threshold)
IRF = IRF/np.sum(IRF)
x = np.asarray([1e-3 * i * res * 2**binning  for i in range(len(IRF))])

np.savetxt('./IRF_th260P.txt',np.stack([x,IRF]),delimiter = ' ')

#plot results
plt.plot(x, IRF)
plt.xlim(0,20)
plt.ylim(1e-4,1e0)
plt.yscale('log')
plt.xlabel('time/ns')
plt.ylabel('counts')
plt.text(2,5e-1, 'FWHM= {} ps'.format(FWHM))
plt.show()

th260.close_all()






