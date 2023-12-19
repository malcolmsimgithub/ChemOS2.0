import pickle
import numpy as np
from tkinter import filedialog
import matplotlib.pyplot as plt
import os, math
from pprint import pprint as pp
import csv

def load_pkl(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data

def calc_freq_spectrum(calc_range = [300,800]): 
    """ 
    calculate the PL spectrum normalized in frequency domain.
    The provided spectrum should be wavelength vs fluence (photons/s/nm)
    The dimension of the calculated spectrum :  second 
    """
    fname = filedialog.askopenfilename()

    data = load_pkl(fname)['PL']['photons']
    l_index = wavelength_to_index(data[0], calc_range[0])

    u_index = wavelength_to_index(data[0], calc_range[1])
    spectrum = data[:, l_index :u_index+1]
    
    area = 0 
    for i, value in enumerate(spectrum[1]):
        if i == 0:
            area += value * (spectrum[0][i+1]-spectrum[0][i])
        elif i == len(spectrum[1])-1:
            area += value * (spectrum[0][i]-spectrum[0][i-1])
        else : 
            area += value * (spectrum[0][i+1] -spectrum[0][i-1])/2
    
    freq_spectrum = 1E9*(spectrum[0]*1E-9)**2*spectrum[1]/(3E8 * area)

    fig =plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(spectrum[0], spectrum[1])
    ax1.set_ylim([-np.max(spectrum[1]) * 0.2, np.max(spectrum[1]) * 1.2])
    ax2 = ax1.twinx()
    ax2.plot(spectrum[0], freq_spectrum, color = 'r')
    ax2.set_ylim([-np.max(freq_spectrum) * 0.2, np.max(freq_spectrum) * 1.2])
    plt.show()

    return np.asarray[spectrum[0], freq_spectrum]

def calc_gain_spectrum(calc_range = [300,800], ksr = 1E9, n = 1.8): 
    """ 
    calculate the normalized PL spectrum in frequency domain.
    The spectrum should be wavelength vs fluence (photons/s/nm)
    The dimension of the calculated spectrum :  second 
    """
    fname = filedialog.askopenfilename()

    data = load_pkl(fname)['PL']['photons']
    l_index = wavelength_to_index(data[0], calc_range[0])
    u_index = wavelength_to_index(data[0], calc_range[1])
    spectrum = data[:, l_index :u_index+1]


    area = 0 
    for i, value in enumerate(spectrum[1]):
        if i == 0:
            area += value * (spectrum[0][i+1]-spectrum[0][i])
        elif i == len(spectrum[1])-1:
            area += value * (spectrum[0][i]-spectrum[0][i-1])
        else : 
            area += value * (spectrum[0][i+1] -spectrum[0][i-1])/2
    
    gain_spectrum = 1E4*1E9*(spectrum[0]*1E-9)**4*ksr*spectrum[1]/(8*math.pi*3E8 * area)

    fig =plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(spectrum[0], spectrum[1])
    ax1.set_ylim([-np.max(spectrum[1]) * 0.2, np.max(spectrum[1]) * 1.2])
    ax2 = ax1.twinx()
    ax2.plot(spectrum[0], gain_spectrum, color = 'r')
    ax2.set_ylim([-np.max(gain_spectrum) * 0.2, np.max(gain_spectrum) * 1.2])
    plt.show()

    return np.asarray([spectrum[0], gain_spectrum])

def wavelength_to_index(data, wavelength):
    index = np.where(data >= wavelength)[0][0]
    return index


def find_abs_end(threshold = 0.01, analysis_range = [300,800]):

    fname = filedialog.askopenfilename()

    data = load_pkl(fname)['absorption']['absorbance']

    l_index, u_index = wavelength_to_index(data[0], analysis_range[0]), wavelength_to_index(data[0], analysis_range[1])

    spectrum = data[:, l_index: u_index+1]
    spectrum[1] = spectrum[1]/np.max(spectrum[1])

    max_index = np.argmax(spectrum[1])
    end_index = None

    for i in range(len(spectrum[1]) - max_index):
        if spectrum[1][max_index + i] < threshold:
            end_index = max_index + i
            break
    if end_index:
        end_wavelength = spectrum[0][end_index]
        print(end_wavelength)
        plt.plot(spectrum[0], spectrum[1])
        plt.show()
        return end_index, end_wavelength
    else:
        raise Exception('absorption spectrum does not go below threshold. Check the baseline')

def res_plot(d_type, plot_range = [300, 800], normalize = False, s_filename=None):

    filelist = filedialog.askopenfilenames()
    for f in filelist:
        data = load_pkl(f)
        d_name = os.path.basename(f).split('.')[0]
        if d_type == 'absorbance':
            d = data['absorption']['absorbance']
            ylabel = 'absorbance'
        elif d_type == 'PL':
            d = data['PL']['energy']
            ylabel = 'intensity/a.u.'
        d_range = [wavelength_to_index(d[0],plot_range[0]), wavelength_to_index(d[0],plot_range[1]+1)]
        if normalize:
            scale = np.max(d[1][d_range[0]:d_range[1]])
        else :
            scale = 1
        plt.plot(d[0][d_range[0]:d_range[1]], d[1][d_range[0]:d_range[1]]/scale, label = d_name)
        plt.ylabel(ylabel)          
        plt.xlabel('wavelength / nm')
        plt.xlim(plot_range)
    plt.legend()
    if s_filename:
        plt.savefig(s_filename)
    plt.show()


def comparison(data_name1, data_name2, plot = True):
    filelist = filedialog.askopenfilenames()
    d = [0 for _ in range(len(filelist))]
    for i, f in enumerate(filelist):
        data = load_pkl(f)
        d[i] = data[data_name1][data_name2]
        print('%s : %s' %(os.path.basename(f), d[i]))

    if plot:
        plt.plot(d)
        plt.show()

def pair_plot(data1, data2):
    filelist = filedialog.askopenfilenames()
    d1 = []
    d2 = []
    for i, f in enumerate(filelist):
        data = load_pkl(f)
        d = data
        for i in range(len(data1)):
            d = d[data1[i]]        
        d1.append(d)
        # de = ''
        # for s in d:
        #     if s.isdigit():
        #         de += s
        # d1.append(float(de))
        d = data
        for i in range(len(data2)):
            d = d[data2[i]]
        d2.append(d)
    plt.plot(d1, d2, linestyle='None', marker="o")
    plt.show()

def stats(analysis_range = [300,800]):

    filelist = filedialog.askopenfilenames()
    for f in filelist:
        data = load_pkl(f)
        absorption = data['absorption']['absorbance']
        PL = data['PL']['energy']
        d_range = [wavelength_to_index(absorption[0],analysis_range[0]), wavelength_to_index(absorption[0],analysis_range[1]+1)]
        abs_max = np.max(absorption[1][d_range[0]:d_range[1]]) 
        lambda_max_abs = absorption[0][np.argmax(absorption[1][d_range[0]:d_range[1]])]
        PL_max = np.max(PL[1][d_range[0]:d_range[1]]) 
        lambda_max_PL = PL[0][np.argmax(PL[1][d_range[0]:d_range[1]])]
        uv_abs = data['uv']['absorbed_power(W)']*1000
        print('max_absorbance : {:.3f} ({:.1f} nm), max_PL : {:.3f} ({:.1f} nm), uv_abs : {:.3f} mW'\
            .format(abs_max, lambda_max_abs, PL_max, lambda_max_PL, uv_abs))
        print('PL_max/abs_max : {:.3f}, PL_max/uv_max : {:.3f}, abs_max/uv_max : {:.3f}'\
            .format(PL_max/abs_max, PL_max/uv_abs, abs_max/uv_abs))
        print('relative_QE : {:.5f}'.format(data['PL']['relative_QE']))


def output_csv(fname):
    filelist = filedialog.askopenfilenames()

    header = [
        'filename',
        'concentration',
        'abs_max',
        'abs_max_wavelength',
        'abs_end',
        'uv_reference', 
        'uv_absorption', 
        'degradation_rate',
        'uv_absorbance_maintenance(at_1min)',
        'PL_max',
        'PL_max_wavelength',
        'relative_QY',
        'max_gain',
        'max_gain_wavelength',
        'PL_lifetime'
    ]

    with open(fname, 'w', newline = '') as f:

        writer = csv.writer(f)
        writer.writerow(header)

        for data_file in filelist:

            d = load_pkl(data_file)
            data = [
                data_file,
                d['metadata']['sample']['concentration(uM)'],
                d['absorption']['max'],
                d['absorption']['lambda_max'],
                d['absorption']['end_wavelength'],
                d['uv']['reference'],
                d['uv']['absorption'],
                d['uv']['degradation_rate'],
                d['uv']['absorption_maintenance(at_1min)'],
                d['PL']['max'],
                d['PL']['lambda_max'],
                d['PL']['relative_QY'],
                d['PL']['max_gain_factor'],
                d['PL']['max_gain_wavelength'],
                d['TE']['fitting_results']['tau']
            ]
            writer.writerow(data)
        



def metadata():

    filelist = filedialog.askopenfilenames()
    for f in filelist:
        data = load_pkl(f)
        pp(data['metadata'])

if __name__=='__main__':
    #res_plot(d_type = 'absorbance', normalize=False)
    #stats()
    #metadata()
    #comparison('uv', 'reference')
    #pair_plot(['metadata', 'sample', 'concentration'],['absorption','max'])
    #pair_plot(['uv', 'reference'],['PL','relative_QY'])
    #pair_plot(['absorption','max'],['PL','max'])
    #pair_plot(['uv','absorption'],['PL','max'])
    #pair_plot(['absorption','max'],['uv','absorption'])
    #calc_freq_spectrum()
    #calc_gain_spectrum()
    #find_abs_end(threshold=0.01)

    fname = 'C:/Users/smily/Dropbox/PythonScript/kazu/data/20200817_test_abs_PL_TE_3/summary3.csv'
    output_csv(fname)