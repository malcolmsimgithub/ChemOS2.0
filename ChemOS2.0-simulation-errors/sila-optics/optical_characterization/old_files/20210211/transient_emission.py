from pylab.instruments import PS5242D, TH260, ThorlabsFW212C
from pylab.instruments.pypowerusb.powerUSB import powerUSB as pUSB
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os, time, copy, csv
from datetime import datetime
import json
from statistics import mean
import pickle

filedir = os.path.dirname(os.path.realpath(__file__))

class Transient_Emission():
    """  
    This is the class to carry out the transient emission experiment by using the following components.
    - TCSPC board : TimeHarp 260 (Picoquant)
    - Oscilloscope : picoscope 5242D (Picotech)
    - USBpowerbar : powerUSB (to switch detector on/off)
    - ND filter wheel : FW212C (Thorlabs)

    Parameters
    --------------------
    filter_visa : str
        visa address of the ND filter wheel
    """
    def __init__(self, device = True, DB = False):

        self._load_setting('%s\.config_transient_emission.dat' %filedir) 
        self.device = device

        if self.device:
            self.fw = ThorlabsFW212C(self.setting['visa_filter'], pos_count=12)
            self.ps = PS5242D('8BIT')
            self.pUSB = pUSB()
            self.th260 = TH260(0)
            self.res = self.th260.get_resolution()
            self._th260_setting()

        if DB:
            from MDB_client import MDB_client_TE
            self.db = MDB_client_TE(config_file = '%s\.config_MDB_client.dat' %filedir, admin=False)


    def _load_setting(self, config_file):
        
        with open(config_file) as content:
            self.config = json.loads(content.read())

        self.setting = {}

        for val in self.config.values():
            for key in val.keys(): 
                if val[key] == 'True':
                    val[key] = True
                elif val[key] == 'False':
                    val[key] = False
                self.setting[key] = val[key]


    def _th260_setting(self):
        self.th260.set_sync_divider(1)
        self.th260.set_cfd(channel='SYNC', cfd_level=self.setting['SYNC_cfd_level'], zero_cross=self.setting['SYNC_zero_cross'])
        self.th260.set_cfd(channel=0,      cfd_level=self.setting['channel_cfd_level'], zero_cross=self.setting['channel_zero_cross'])

        self.th260.set_channel_offset(channel='SYNC', offset=self.setting['SYNC_offset'])
        self.th260.set_channel_offset(channel=0,      offset=self.setting['channel_offset'])

        self.th260.set_binning(self.setting['binning'])
        self.th260.set_offset(self.setting['offset'])


    def laser_on(self, frequency):   #frequency in Hz, amplitude in mV
        """  
        Turn on laser by using picoscope as a external trigger source.
        The laser driver detects the rising edge of the trigger signal.
        The driver should be in external triggering mode to use this function. 

        Parameters
        --------------------
        frequency : float
            laser frequency in Hz
        """
        self.ps.setBuiltInSignal(1, frequency, 1500, 3000)
        self.ps.startSignal()


    def laser_off(self):
        """  
        Turn off laser by using picoscope as a external trigger source.
        The driver should be in external triggering mode to use this function. 
        """
        self.ps.stopSignal()


    def detector_on(self, port = 1):
        """  
        Turn on detector by using the power usb.

        Parameters
        --------------------
        port : int
            port number where the detector connected
        """
        self.pUSB.setport(port, 1)
        # self.pUSB.setport_all(1, 0, 0)
        while True:
            if self.th260.get_count_rate(0) < 1:
                time.sleep(1)
            else: 
                print('detector:on')
                break


    def detector_off(self, port = 1):
        """  
        Turn off detector by using the power usb.

        Parameters
        --------------------
        port : int
            port number where the detector connected
        """
        self.pUSB.setport_all(0, 0, 0)
        print('detector:off')


    def set_filter(self, position): 
        """  
        set position of the ND filter

        Parameters
        --------------------
        position : int
            filter position(0-12)
        """     
        self.fw.set_position(position)


    def adjust_filter(self, min_rate, max_rate, s_position = 12): 
        """  
        Adjust the ND filter based on the detected light intensity.
        The filter will be set to the position where the itensity get
        closest to the max_rate (< max_rate).
        If the intensity is less than min_rate at the no filter position (position 1)
        or more than max_rate at the strongest filter position(position 12),
        the function raises the exception.

        Parameters
        --------------------
        min_rate : float
            minimum count rate (relative to the sync rate)
        max_rate : float
            maximum count rate (relative to the sync rate)
        s_position : int
            starting position of the filter (1-12)
        """  
        flg = 0
        self.set_filter(s_position)
        self.laser_on(1e6)
        time.sleep(1)
        while True:
            pos = self.fw.get_position()
            s_count, c_count = self.th260.get_sync_rate(), self.th260.get_count_rate(0)
            print(pos, s_count, c_count)
            if c_count <  max_rate * s_count:
                if pos == 1:
                    if c_count < min_rate * s_count:
                        print('Warning : PL intensity did not reach minimum count rate')
                        flg  = 1
                        break
                    else:
                        break
                else:
                    self.set_filter(pos-1)
            else: 
                if pos == 12:
                    print('Warning : Falsed to reduce PL intensity to less than maximum count rate')
                    flg = 1
                    break
                else:
                    self.set_filter(pos+1)
                    break
        self.laser_off()
        return flg

    def measure(self, overflow, frequency, timeout=60):
        """  
        Start accumulating the signal.

        Parameters
        --------------------
        overflow : int
            accumulation stop when the maximum count reaches this value
        frequency : int
            frequency of the laser in Hz.
            The laser driver should be in the external triggering mode 
        timeout : int
            the maximum accumulation time in seconds
        
        Returns
        --------------------
        x : 1D numpy array
            time in nanoseconds
        hist : 1D numpy array
            detected count
        """ 
        self.th260.set_overflow(overflow)
        self.laser_on(frequency)
        time.sleep(1)
        print(self.th260.get_sync_rate())
        hist = np.asarray(self.th260.measure(timeout))
        self.laser_off()
        x = np.asarray([1e-3 * i * self.res * 2**self.setting['binning']  for i in range(len(hist))])
        return x, hist


    def exp_fit(self, x, hist, weight = 1, rise = 50):
        """  
        Do exponential fitting to the measurement results.

        Parameters
        --------------------
        x : 1D numpy array
            time in nanoseconds
        hist : 1D numpy array
            detected count
        weight : float or str
            float :adjust the weight of the data points by; 
                (amp * np.exp(-1*(x-x0)/tau))**weight
            'sqrt' : weight by square root of y
            'lin' : weight by y
            'inv' : weight by 1/y
        rise : int
            number of data point for rising edge used to calculate the background.
            (bg = np.average(hist[0:max_index-rise]))

        Returns
        --------------------
        fit_result : dict
            {'amp', 'tau', 'x0', 'bg', 'max_index' , 'covar'}
        """ 

        def exp(x, amp, tau, x0):
            return (amp * np.exp(-1*(x-x0)/tau))**wt

        if type(weight) == int or type(weight) == float:
            wt = weight
            sigma = None
        elif type(weight) == str:
            wt = 1

        c_max, max_index = np.max(hist), np.argmax(hist)
        bg = np.average(hist[0:max_index-rise])
        a_hist = hist[max_index:]-bg
        a_hist[a_hist<0]=0.00001   
        a_x, a_hist = x[max_index:]-x[max_index], (a_hist)**wt
        init_vals = [c_max, 10, 0]

        if weight == 'sqrt':
            sigma = a_hist**0.5
        elif weight == 'lin':
            sigma = a_hist
        elif weight == 'inv':
            sigma = a_hist**(-1)
        else:
            sigma = None

        fit_vals, covar = curve_fit(exp, a_x, a_hist, sigma = sigma, \
                             absolute_sigma = True, p0=init_vals, bounds=(0,[np.inf, np.inf, 4]))

        rss = np.sum((a_hist - exp(a_x, fit_vals[0], fit_vals[1], fit_vals[2]))**2)
        tss = np.sum((a_hist-np.mean(a_hist))**2)
        r_square = 1 - rss/tss

        print(fit_vals, bg, r_square)
        fit_results ={'amp': fit_vals[0], 'tau' : fit_vals[1], 'x0' : fit_vals[2], 
                                                'bg' : bg, 'max_index' : int(max_index), 'R2': r_square}
        
        return fit_results



    def exp_fit_log(self,x, hist, rise = 50):
        """  
        Do exponential fitting to the logalism of the measurement results 

        Parameters
        --------------------
        x : 1D numpy array
            time in nanoseconds
        hist : 1D numpy array
            detected count
        rise : int
            number of data point for rising edge used to calculate the background.
            (bg = np.average(hist[0:max_index-rise]))

        Returns
        --------------------
        fit_result : dict
            {'amp', 'tau', 'x0', 'bg', 'max_index' , 'covar'}
        """ 
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

        rss = np.sum((a_hist - exp_log(a_x, fit_vals[0], fit_vals[1], fit_vals[2]))**2)
        tss = np.sum((a_hist-np.mean(a_hist))**2)
        r_square = 1 - rss/tss

        print(fit_vals, bg, r_square)
        fit_results ={'amp': fit_vals[0], 'tau' : fit_vals[1], 'x0' : fit_vals[2], 
                                                'bg' : bg, 'max_index' : int(max_index), 'R2': r_square}

        return fit_results


    def double_exp_fit(self, x, hist, weight = 1, rise = 50):
        """  
        Do double exponential fitting to the measurement results.

        Parameters
        --------------------
        x : 1D numpy array
            time in nanoseconds
        hist : 1D numpy array
            detected count
        weight : float or str
            float :adjust the weight of the data points by; 
                (amp * np.exp(-1*(x-x0)/tau))**weight
            'sqrt' : weight by square root of y
            'lin' : weight by y
            'inv' : weight by 1/y
        rise : int
            number of data point for rising edge used to calculate the background.
            (bg = np.average(hist[0:max_index-rise]))

        Returns
        --------------------
        fit_result : dict
            {'amp1', 'tau1', 'amp2', 'tau2', x0', 'bg', 'max_index' , 'covar'}
        """ 

        def double_exp(x, amp1, amp2, tau1, tau2, x0):
            return (amp1 * np.exp(-1*(x-x0)/tau1) + amp2 * np.exp(-1*(x-x0)/tau2))**wt

        if type(weight) == int or type(weight) == float:
            wt = weight
            sigma = None
        elif type(weight) == str:
            wt = 1

        if weight == 'sqrt':
            sigma = a_hist**0.5
        elif weight == 'lin':
            sigma = a_hist
        elif weight == 'inv':
            sigma = a_hist**(-1)
        else:
            sigma = None

        c_max, max_index = np.max(hist), np.argmax(hist)
        bg = np.average(hist[0:max_index-rise])
        a_hist = hist[max_index:]-bg
        a_hist[a_hist<0]=0.00001   
        a_x, a_hist = x[max_index:]-x[max_index], a_hist**wt
        init_vals = [c_max, c_max/5, 5, 20, 0]

        fit_vals, covar = curve_fit(double_exp, a_x, a_hist, sigma = sigma, \
                             absolute_sigma = True, p0=init_vals, bounds=(0,[np.inf, np.inf, np.inf, np.inf, 4]))

        rss = np.sum((a_hist - double_exp(a_x, fit_vals[0], fit_vals[1], fit_vals[2], fit_vals[3], fit_vals[4]))**2)
        tss = np.sum((a_hist-np.mean(a_hist))**2)
        r_square = 1 - rss/tss

        print(fit_vals, bg, r_square)
        fit_results ={'amp1': fit_vals[0], 'tau1' : fit_vals[2], 'amp2': fit_vals[1], 'tau2' : fit_vals[3], 'x0' : fit_vals[4], 
                                                'bg' : bg, 'max_index' : int(max_index), 'R2': r_square}

        return fit_results


    def plot_raw(self, x, hist, freq):
        """  
        plot raw result

        Parameters
        --------------------
        x : 1D numpy array
            time in nanoseconds
        hist : 1D numpy array
            detected count
        freq : int
            laser frequency for the measurement
        """ 
        plt.plot(x, hist)
        plt.xlim(10,0.9*1e9/freq)
        plt.ylim(1e0,2e4)
        plt.yscale('log')
        plt.xlabel('time/ns')
        plt.ylabel('counts')
        plt.show()


    def plot_fit(self, data, fname = None, plot = False):
        """  
        plot raw result with fitting results.

        Parameters
        --------------------
        x : 1D numpy array
            time in nanoseconds
        data : dict
            result dictionary
        freq : int
            laser frequency for the measurement
        fname : str, or None
            save plot if not None(.png)  
        plot : bool
            show plot if True
        """ 
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, xlim = [10,0.9*1e9/data['metadata']['excitation_frequency']],\
                ylim = [1e0, np.max(data['raw_data'][1])*2], yscale = 'log', xlabel = 'time/ns', ylabel = 'counts')
        ax.plot(data['raw_data'][0], data['raw_data'][1])
        if 'fitting_results' in data.keys():
            if 'tau' in data['fitting_results'].keys():
                ax.plot(data['fitted_data'][0], data['fitted_data'][1])
                ax.text(0.7, 0.85, 'tau = {:.2f} ns'.format(data['fitting_results']['tau']), transform=ax.transAxes)
                ax.text(0.7, 0.77, 'R2 = {:.4f}'.format(data['fitting_results']['R2']), transform=ax.transAxes)
            if 'tau1' in data['fitting_results'].keys():
                ax.plot(data['fitted_data'][0], data['fitted_data'][1])
                ax.text(0.7, 0.9, 'tau1 = {:.2f} ns'.format(data['fitting_results']['tau1']), transform=ax.transAxes)
                ax.text(0.7, 0.825, 'amp1 = {:.4f}'.format(data['fitting_results']['amp1']), transform=ax.transAxes)
                ax.text(0.7, 0.75, 'tau2 = {:.2f} ns'.format(data['fitting_results']['tau2']), transform=ax.transAxes)
                ax.text(0.7, 0.675, 'amp2 = {:.4f}'.format(data['fitting_results']['amp2']), transform=ax.transAxes)
                ax.text(0.7, 0.6, 'R2 = {:.4f}'.format(data['fitting_results']['R2']), transform=ax.transAxes)
        #plt.title('%s uM %s in %s' %(data['metadata']['sample']['concentration(uM)'], \
        #                        data['metadata']['sample']['name'], data['metadata']['sample']['solvent']))
        if fname:
            print(fname + "_TE.png")
            plt.savefig(fname + '_TE.png')
        if plot:
            plt.show()

        plt.close()


    def save_pkl(self, filename, data):
        """  
        Save measurement results in pkl file.
        {'raw' : raw_data(2D numpy array), 
         'fit' : fittting_result(2D numpy array), 
         'fit_param' : fitting parameters(dict)}

        Parameters
        --------------------
        filename : str
            filename to be saved
        data : dict
            result dictionary
        """ 
        with open(filename, 'wb') as f:
            pickle.dump(data, f)


    def load_pkl(self, filename):
        """  
        load pkl file.

        Parameters
        --------------------
        filename : str
            filename to be saved

        Returns
        --------------------
        data : dict
        """ 
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data


    def close_all(self):
        """  
        Close timeharp, picoscope and detector
        """ 
        self.th260.close_all()
        self.ps.close()
        #self.detector_off()


    def _transpose_list(self, l):

        max_dnum = max([len(row) for row in l])

        t_list = []
        for i in range(max_dnum):
            t_list.append([])
            for row in l:
                if i <= len(row) -1:
                    t_list[i].append(row[i])
                else: t_list[i].append(None)
        
        return t_list


    def _create_result_dict(self, data, setting, freq, save_filename):
        
        config = copy.deepcopy(self.config)

        for key, val in config.items():
            for k in val.keys():
                config[key][k] = setting[k]

        config['experiment_timestamp'] =str(datetime.now())
        config['excitation_frequency'] = freq
        if save_filename:
            config['filename'] = '%s.pkl' %save_filename
        else: 
            config['filename'] = None

        result = {'raw_data' : data['raw_data'],
                  'metadata' : config}

        if 'fitted_data' in data.keys():
            result['fitted_data'] = data['fitted_data']
            result['fitting_results'] = data['fitting_results']

        return result


    def do_fitting(self, x, hist, fit_order = 'auto', fit_weight = 1):
        """  
        Do fitting on the data. Fitting is done either single-exponential or double-exponential 
        depends on the fit_order setting.
        If the fit_order is set as auto, single or double which gives the smaller R2 value will be selected.
        
        Parameters
        --------------------
        x : list
        hist : list 
            histgram data
        fit_order :int or str
            fitting order (1 or 2 or auto).
        fit_weight : float or str 
            float :adjust the weight of the data points by; 
                (amp * np.exp(-1*(x-x0)/tau))**weights
            str : 'sqrt', 'lin', inv'

        Returns
        --------------------
        data: dict
            {'raw_data' : 2D numpy array, 
             'fitted_data' : 2D numpy array, 
             'fitting_results' : dict}
        """ 
        fit1_R2, fit2_R2 = 0, 0

        if fit_order== 1 or fit_order== 'auto':
            try:
                fit1_res = self.exp_fit(x, hist, weight = fit_weight) 
                fit1_R2 = fit1_res['R2']   
            except:
                print('Fitting error (1st order)')
        if fit_order== 2 or fit_order== 'auto':
            try:
                fit2_res = self.double_exp_fit(x, hist, weight = fit_weight)
                fit2_R2 = fit2_res['R2']   
            except:
                print('Fitting error (2nd order)')

        if fit1_R2 == 0 and fit2_R2 == 0:
            print('Fitting error! Only the raw data will be saved')
            result = {'raw_data': np.vstack((x, hist))}

        elif fit1_R2 >= fit2_R2 - 0.0002:       
            x_fit = x[fit1_res['max_index']:]
            hist_fit = [fit1_res['amp']*np.exp(-1*(xx-x[fit1_res['max_index']]-fit1_res['x0'])/fit1_res['tau'])+fit1_res['bg'] for xx in x_fit]

            result = {'raw_data': np.vstack((x, hist)), 'fitted_data': np.vstack((x_fit, hist_fit)), 'fitting_results' : fit1_res}

        else:
            x_fit = x[fit2_res['max_index']:]
            hist_fit = [fit2_res['amp1']*np.exp(-1*(xx-x[fit2_res['max_index']]-fit2_res['x0'])/fit2_res['tau1'])+\
                fit2_res['amp2']*np.exp(-1*(xx-x[fit2_res['max_index']]-fit2_res['x0'])/fit2_res['tau2'])+fit2_res['bg'] for xx in x_fit]

            result = {'raw_data': np.vstack((x, hist)), 'fitted_data': np.vstack((x_fit, hist_fit)), 'fitting_results' : fit2_res}
        
        return result



    def measure_TE(self, save_filename = None, **kwargs):
        """  
        Carry out the transient emission measurement based on the setting parameters.
        Fitting is done either single-exponential or double-exponential depends on the fit_order setting.
        If the fit_order is set as auto, single or double which gives the smaller R2 value will be selected. 

        Parameters
        --------------------
        save_filename : str, or None 
            save filename for the results without extention
        min_rate : float 
            minimum count rate necessary for the measurement)
        max_rate : float 
            maximum count rate for the measurement
        accumeration : int 
            maximum count to be accumulated
        fit_order :int or str
            fitting order (1 or 2 or auto).
        fit_weight : float or str 
            float :adjust the weight of the data points by; 
                (amp * np.exp(-1*(x-x0)/tau))**weights
            str : 'sqrt', 'lin', inv'
        initial_filter_position' : int 
        initial_freq :  int 
            initial laser frequency
        do_plot : bool
            show result plot if True
        name : str
            sample name
        concentration : str
            concentration of the sample
        solvent : str
            solvent used

        Returns
        --------------------
        results: dict
            {'raw_data' : 2D numpy array, 
             'fitted_data' : 2D numpy array, 
             'fitting_results' : dict,
             'metadata' : dict}
        """ 
        #TODO make dir

        setting = self.setting

        for key in kwargs.keys():
            if key in setting.keys():
                setting[key] = kwargs[key]

        flg = self.adjust_filter(setting['min_rate'], setting['max_rate'], s_position = setting['initial_filter_position']) 
        if flg == 1:
            self.laser_off()
            return

        # set_frequency
        x, hist = self.measure(100, setting['initial_frequency'])
        fit_res = self.exp_fit(x, hist, weight = setting['fit_weight'])
        freq = int(min(1/(fit_res['tau']*9.212*3*1e-9), 2e7))
        print('measurement is carried out at %s Hz' %freq)
        
        #Measure Histgram 
        x, hist = self.measure(setting['accumulation'], frequency = freq)
    
        #fitting
        fit_res = self.do_fitting(x, hist, fit_order = setting['fit_order'], fit_weight = setting['fit_weight'])
       
        results = self._create_result_dict(fit_res, setting, freq, save_filename)

        #plot results
        self.plot_fit(results, fname = save_filename, plot = setting['do_plot'])        
        
        #self.close_all()

        #export results
        if save_filename:
            self.save_pkl(save_filename + '_TE.pkl', results)

            l = results['raw_data'].tolist()
            if 'fitting_results' in results.keys():
                l.append(results['fitted_data'].tolist())

            with open(save_filename + '_TE.csv', 'w', newline = '') as f:
                writer = csv.writer(f)
                for row in self._transpose_list(l):
                    writer.writerow(row)

            # np.savetxt(save_filename + '_TE.csv', np.transpose(results['raw_data']))
            # if 'fitting_results' in results.keys():
            #     np.savetxt(save_filename + '_fit.csv', np.transpose(results['fitted_data']))

        return results
  
  
    def measure_TE_emu(self, data, freq, save_filename = None, **kwargs):
        """  
        emurator for fitting test
        """

        setting = self.setting

        for key in kwargs.keys():
            setting[key] = kwargs[key]

        x, hist = data[0], data[1]
    
        #fitting
        fit_res = self.do_fitting(x, hist, fit_order = setting['fit_order'], fit_weight = setting['fit_weight'])
       
        results = self._create_result_dict(fit_res, setting, freq, save_filename)

        #plot results
        self.plot_fit(results, fname = save_filename, plot = setting['do_plot'])        

        #export results
        if save_filename:
            self.save_pkl(save_filename + '.pkl', results)

            l = results['raw_data'].tolist()
            if 'fitting_results' in results.keys():
                l.append(results['fitted_data'].tolist())

            with open(save_filename + '_TE.csv', 'w', newline = '') as f:
                writer = csv.writer(f)
                for row in _transpose_list(l):
                    writer.writerow(row)
            # np.savetxt(save_filename + '.csv', np.transpose(np.vstack((x, hist))))
            # if 'fitting_results' in data.keys():
            #     np.savetxt(save_filename + '_fit.csv', np.transpose(np.vstack((x_fit, hist_fit))))
        
        return results
  
if __name__ == "__main__":

    

    ###measurement##########
    """
    experimental parameters are in ".config_transient_emission.dat" 
    """

    # TE = Transient_Emission(device = True, DB = False)
    # TE.detector_on()

    # sample_info ={
    #     'name' : 'di-nBu-PTPTP',
    #     'concentration' : '2uM',
    #     'solvent' : 'THF'
    # }

    # path = 'C:/Users/MatterLab/Dropbox/PythonScript/kazu/data/20200714_test_autoPL/'
    # sample_name = '20uL-di-nBu-PTPTP_in_THF' 

    # #sfilename = '%s/test.pkl' %filedir
    # sfilename = None
    # result = TE.measure_TE(save_filename = path + sample_name, fit_order = 'auto', fit_weight = 0.7, do_plot = True, **sample_info)
  
    ###emulator##########
  
    TE = Transient_Emission(device = False, DB = False)

    fname = 'C:/Users/smily/Dropbox/PythonScript/kazu/data/for_automation/sample_measured/mix1_100uM_char_diphenylanthracene_9.313_diphenylanthracene_PL/mix1_100uM_char_diphenylanthracene_9.313_diphenylanthracene_PL.pkl'
    
    data = TE.load_pkl(fname)
    print(data.keys())

    TE.measure_TE_emu(data['raw_data'], 7000000, save_filename = None, fit_order = 1, fit_weight = 0.7, do_plot = True) 
  
  
  
  
 
