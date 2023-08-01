import pylab.instruments as ins
from pylab.manager import Manager, Logger
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy import signal
import pickle
import math
import copy
import threading

default_LEDs = {
    '365' : {'driver' : 1, 'initial_power' : 10},
    '415' : {'driver' : 2,  'initial_power' : 10}
}

class Abs_PL():
    
    def __init__(self, leds = None, correction_datafile = None, device = True):
        
        self.device = device

        if self.device:

            self.led_list = default_LEDs
            if leds:
                for key, val in self.led_list.items():
                    if key in leds.keys():
                        val['initial_power'] = leds[key]
                    else : 
                        val['initial_power'] = 0
            print(self.led_list)
            self.led_driver = ins.ThorlabsDC4100("ASRL9::INSTR")
            self.leds = {key : self.led_driver.led[val['driver']] for key, val in LEDs.items()}
            # self.led_365 = self.led_driver.led[1]
            self.shutter = ins.ThorlabsKSC101("68001529")

            self.lamp = ins.DH_mini('ASRL10::INSTR', verbose=False)
            self.spectrometer = ins.ThorlabsCCS("USB0::0x1313::0x8089::M00562925::RAW")
            self.powermeter = ins.ThorlabsPM100D("USB0::0x1313::0x8078::P0025495::INSTR", averages=2000)

            self.wl, _ = self.spectrometer.get_wavelength_data()
            self.dark_spectrum = []
            
            for key, val in self.led_list.items():
                val['wavelength'] = self.leds[key].get_wavelength()
                val['led_photon_energy'] = 2.998E8 * 6.626E-34 /(val['wavelength']*1E-9)  #photon energy in J
            # self.led_wl = self.led_365.get_wavelength()
            # self.led_photon_energy = 2.998E8 * 6.626E-34 /(self.led_wl*1E-9)  #photon energy in J

            #TODO correction data
            if correction_datafile:
                self.correction_data = self._load_correction_data(correction_datafile)
            else:
                self.correction_data = np.ones(len(self.wl))

            self._initialize_device()


    def _initialize_device(self):

        print('***initializing the device***********')
        self.lamp_on()
        self.lamp_shutter_close()

        for key, val in self.led_list.items():
            if val['initial_power'] != 0:
                self.led_on(key, val['initial_power'])
            else:
                self.led_off(key)
        self.led_shutter_close()
        ###TODO check
        self.powermeter.set_mode('POWER')
        self.powermeter.set_autorang('ON')  #activate powermeter auto-range

        print('***initilizing done***********')


    def _load_correction_data(self, filename):
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        return data

    ###white light functions#####
    def lamp_on(self):
        self.lamp.deuterium_on()
        self.lamp.halogen_on()

    def lamp_off(self):
        self.lamp.deuterium_off()
        self.lamp.halogen_off()

    def lamp_shutter_open(self):
        self.lamp.shutter_open()

    def lamp_shutter_close(self):
        self.lamp.shutter_close()


    def check_lamp_status(self):

        status = self.lamp.check_status()

        for l in ['halogen', 'deuterium']:
            if status[l] == 'off':
                print('causion : %s was off. The lamp might be not stable yet' %l)
        
        return status


    ###led functions#############
    def led_on(self, wavelength, percentage, off_others = False):
        if off_others:
            for key in self.leds.keys():
                if key != str(wavelength):
                    self.led_off(key)
        self.leds[str(wavelength)].set_percentage(percentage)
        self.leds[str(wavelength)].on()

    def led_off(self, wavelength):
        self.leds[str(wavelength)].off()

    def led_shutter_open(self):
        self.shutter.open()
    
    def led_shutter_close(self):
        self.shutter.close()

    def check_led_state(self, wavelength):
        
        status = self.leds[str(wavelength)].get_state()

        if status == 'OFF':
            print('causion : led <%s> was off. The led might be not stable yet' %wavelength)

        return status

    ###power meter functions#####
    def measure_powers(self, counts):

        powers = []
        for _ in range(counts):
            power = self.powermeter.measure()
            powers.append([time.time(), power])

        return powers


    def average_power(self, counts):
        power = 0.
        for _ in range(counts):
            power += self.powermeter.measure()
        return power/counts


    def measure_uv_power(self, wavelength, counts, led_power = 80):

        self.check_led_state(wavelength)
        self.led_on(wavelength, percentage = led_power, off_others=True)

        self.lamp_shutter_close()
        self.led_shutter_open()
        time.sleep(1)
        power = self.average_power(counts)
        self.led_shutter_close()

        return power

    ###spectrometer functions#####
    def measure_uv_absorption(self, wavelength, reference, counts, led_power = 80):

        power = self.measure_uv_power(wavelength, counts, led_power)

        absorbed_power = reference - power
        absorbed_photon = absorbed_power/self.led_list[str(wavelength)]['led_photon_energy']  #photons/s

        return absorbed_power, absorbed_photon 



    def measure_spectrum(self, repeats, integration_time, filter_size = None, dark_correction = True, spectral_correction = False, do_plot = False):
        
        spec, _ = self.spectrometer.measure(repeats = repeats, integration_time = integration_time)
        
        if filter_size:
            fil = np.ones(filter_size)/filter_size
            spec = signal.convolve(spec, fil, mode='same')

        if dark_correction:
            if len(self.dark_spectrum) == 0:
                raise Exception('No dark spectrum stored. Measure it by "measure_dark_spectrum".')
            spec = spec - self.dark_spectrum

        if spectral_correction:
            spec = spec / self.correction_data

        if do_plot:
            plt.plot(self.wl, spec)
            plt.xlabel('wavelength / nm')
            plt.ylabel('intensity / a.u.')
            plt.show()
        
        return np.asarray([self.wl, spec])


    def measure_spectrum_trace(self, repeats, integration_time, filter_size = None, dark_correction = True, spectral_correction = False, do_plot = False):
        
        results, _ = self.spectrometer.measure_trace(repeats = repeats, integration_time = integration_time)

        if filter_size:
            fil = np.ones(filter_size)/filter_size
            results['average'] = signal.convolve(results['average'], fil, mode='same')
            for i, trace in enumerate(results['time_trace']):
                results['time_trace'][i] = signal.convolve(trace, fil, mode='same')

        if dark_correction:
            if len(self.dark_spectrum) == 0:
                raise Exception('No dark spectrum stored. Measure it by "measure_dark_spectrum".')
            results['average'] =results['average'] - self.dark_spectrum
            for i, trace in enumerate(results['time_trace']):
                results['time_trace'][i] = trace - self.dark_spectrum

        if spectral_correction:
            results['average'] =results['average'] / self.correction_data
            for i, trace in enumerate(results['time_trace']):
                results['time_trace'][i] = trace / self.dark_spectrum

        if do_plot:
            plt.plot(self.wl, results['average'])
            for trace in results['time_trace']:
                plt.plot(self.wl, trace)
            plt.xlabel('wavelength / nm')
            plt.ylabel('intensity / a.u.')
            plt.show()

        results['wavelength'] = self.wl

        return results


    def measure_dark_spectrum(self, repeats, integration_time, filter_size = None, do_plot = False):

        self.lamp_shutter_close()
        self.led_shutter_close()
        
        time.sleep(1)
        spec, _ = self.spectrometer.measure(repeats = repeats, integration_time = integration_time)

        if filter_size:
            fil = np.ones(filter_size)/filter_size
            spec = signal.convolve(spec, fil, mode='same')

        if do_plot:
            plt.plot(self.wl, spec)
            plt.xlabel('wavelength / nm')
            plt.ylabel('intensity / a.u.')
            plt.show()

        self.dark_spectrum = spec



    def measure_transmission_spectrum(self,  repeats, integration_time, filter_size = None, dark_correction = True, do_plot = False):

    
        self.lamp.shutter_open()
        self.lamp_on()

        #if 0 in self.check_lamp_status():        
            #self.lamp_on()

        self.led_shutter_close()

        time.sleep(1)

        transmission_spectrum = self.measure_spectrum(repeats, integration_time, filter_size = filter_size, dark_correction = dark_correction, do_plot = do_plot)

        self.lamp.shutter_close()

        return transmission_spectrum


    # def equilibration(self, repeats, integration_time, equilibration_time, threshold, wavelength, time_step, time_out, filter_size = None):

    #     self.lamp.shutter_open()
    #     self.lamp_on()

    #     self.led_shutter_close()

    #     time.sleep(1)

    #     start_time = time.time()
    #     index = np.where(self.wl >= wavelength)[0][0]
    #     #print(index)

    #     plt.ion()
    #     time_list = []
    #     t = []
    #     i = 0
    #     t_0 = self.measure_spectrum(repeats, integration_time, filter_size = filter_size, dark_correction = False, do_plot = False)[1][index]

    #     time_list.append(0)
    #     t.append(1)

    #     time_ref, t_ref = time_list[0], t[0]

    #     time.sleep(time_step)

    #     while True:

    #         i += 1
            
    #         time_list.append(time.time() - start_time)
    #         t.append(self.measure_spectrum(repeats, integration_time, filter_size = filter_size, dark_correction = False, do_plot = False)[1][index]/t_0)

    #         plt.plot(time_list[0:i], t[0:i])
    #         plt.draw()
    #         plt.pause(0.01)
    #         plt.clf()

    #         print('t : %s, time : %s, e_time : %s' %(t[i], time_list[i], time_list[i]-time_ref))

    #         if time_list[i] > time_out:
    #             print("time out before equilibration")
    #             break

    #         if abs(t[i]-t_ref) < threshold:
    #             if time_list[i] - time_ref > equilibration_time:
    #                 print("equilibration reached")
    #                 break
    #         else:
    #             time_ref, t_ref = time_list[i], t[i]
    #             print('t_ref : %s, time_ref : %s' %(t_ref, time_ref))

    #         time.sleep(time_step)

    #     plt.close()
    #     plt.ioff()
    #     self.lamp.shutter_close()


    def measure_absorption_spectrum(self,  repeats, integration_time, reference, filter_size = None, dark_correction = True, do_plot = False):

        transmission_spectrum = self.measure_transmission_spectrum(repeats, integration_time, filter_size = filter_size, dark_correction = dark_correction, do_plot = False)

        transmittance = transmission_spectrum[1]/reference[1]
        absorbance = - np.log10(transmittance) 

        if do_plot:
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
            ax1.plot(self.wl, absorbance)
            ax2 = ax1.twinx()
            ax2.plot(self.wl, transmittance, color = 'r')
            plt.xlim(300,800)
            plt.show()

        return np.asarray(transmission_spectrum),  np.asarray([self.wl, transmittance]), np.asarray([self.wl, absorbance])


    def adjust_PL_exposure(self, wavelength, initial_exposure, max_exposure, target_intensity, led_power, filter_size = 20, average = 5, max_iter = 5):

        saturation = True
        exposure = initial_exposure

        i = 0

        while saturation and i < max_iter:

            i += 1

            PL_spectrum = self.measure_PL_spectrum(wavelength, average, exposure, led_power=led_power, filter_size=filter_size, dark_correction=False)

            peak_intensity = np.max(PL_spectrum[1])

            if peak_intensity >= target_intensity and peak_intensity < 0.9:
                saturation = False
            elif peak_intensity < target_intensity:
                saturation = False
                exposure =  min(exposure * target_intensity/peak_intensity, max_exposure)
            elif peak_intensity >= 0.9 and peak_intensity < 1:
                saturation = False
                exposure *= target_intensity/peak_intensity
            else:
                exposure *= 0.1
            print('peak intensity : %s, exposure :%s s' %(peak_intensity, exposure))

        return exposure


    def measure_PL_spectrum(self, wavelength, repeats, integration_time, led_power = 80, filter_size = None, dark_correction = True, spectral_correction = True, do_plot = False):

        # if self.check_led_state(wavelength) == 'OFF':
        self.check_led_state(wavelength)
        
        self.led_on(wavelength, led_power, off_others=True)
    
        self.lamp_shutter_close()
        self.led_shutter_open()

        time.sleep(1)

        PL_spectrum = self.measure_spectrum(repeats, integration_time, filter_size = filter_size, \
                             dark_correction = dark_correction, spectral_correction=spectral_correction, do_plot = do_plot)

        return PL_spectrum


    def measure_PL_uv(self, wavelength, repeats, integration_time, led_power = 80, filter_size = None, dark_correction = True, spectral_correction = True, do_plot = False):

        def measure_PL(results, flg):
            PL_results = self.measure_spectrum_trace(repeats, integration_time, filter_size = filter_size,\
                             dark_correction = dark_correction, spectral_correction=spectral_correction, do_plot = do_plot)
            flg['flg'] = False
            results['PL'] = PL_results #{'time', 'average', 'wavelength', 'time_trace'}
            print('PL measurement is done')
        
        def measure_uv(results_uv, flg):
            power = []
            m_time= []
            counts = 0
            while flg['flg'] == True:
                m_time.append(time.time())
                power.append(self.powermeter.measure())
                counts += 1
            results_uv['uv_average'] =  sum(power)/counts   #W
            results_uv['duration'] = m_time[-1] - m_time[0]
            results_uv['uv_time_trace'] = np.asarray([m_time, power])
            results_uv['wavelength'] = wavelength

            print('uv measurement is done')

        # if self.check_led_state() == 'OFF':
        self.check_led_state(wavelength)
        self.led_on(wavelength, led_power, off_others=True)
    
        self.lamp_shutter_close()
        self.led_shutter_open()

        time.sleep(3)

        results, results_uv = {}, {}
        flg = {'flg' : True}

        thread1 = threading.Thread(target=measure_PL, args=([results, flg]))
        thread2 = threading.Thread(target=measure_uv, args=([results_uv, flg]))
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()
        
        results.update(results_uv)
        
        self.led_shutter_close()

        return results


    # def PL_uv_analysis(self, PL_uv_results, uv_reference, PL_exposure, calc_range):

    #     l_index, u_index = self._to_index_range(self.wl, calc_range)
    #     uv, PL = {}, {}

    #     uv['reference'] = uv_reference
    #     uv['absorption'] = uv['reference'] - PL_uv_results['uv_average']
    #     uv['absorbance_time_trace'] = np.asarray([PL_uv_results['uv_time_trace'][0]-PL_uv_results['uv_time_trace'][0, 0], -np.log10(PL_uv_results['uv_time_trace'][1]/uv['reference'])])
    #     PL['energy'] = np.asarray([PL_uv_results['PL']['wavelength'], PL_uv_results['PL']['average']/PL_exposure])
    #     PL['photons'] = np.asarray([PL_uv_results['PL']['wavelength'], PL['energy'][1]*PL['energy'][0]*1E-9])
    #     PL['freq_spectrum'] = self.calc_freq_spectrum(PL['photons'], calc_range= calc_range)
    #     PL['gain_spectrum'] = self.calc_gain_spectrum(PL['photons'], calc_range= calc_range)
    #     PL['relative_QY'] = np.sum(PL['photons'][1][l_index:u_index+1])/uv['absorption']
    #     PL['max'], PL['lambda_max'] = self.find_max(PL['energy'], analysis_range=calc_range)
    #     max_index = np.argmax(PL['energy'][1][l_index:u_index+1])+l_index
    #     PL['time_trace'] = PL_uv_results['PL']['time_trace']
    #     PL['peak_time_trace'] = np.asarray([PL_uv_results['PL']['time'], [np.mean(trace[max_index-3:max_index+3]) for trace in PL['time_trace']]])
        
    #     print('PL_max :  {:.4f} ({:.1f} nm)'.format(PL['max'], PL['lambda_max']))
    #     print('uv_ref : {:.3f} mW'.format(uv['reference'] * 1000))
    #     print('uv_start : {:.3f} mW'.format(PL_uv_results['uv_time_trace'][1, 0] * 1000))
    #     print('uv_end : {:.3f} mW'.format(PL_uv_results['uv_time_trace'][1, -1] * 1000))
    #     print('relative_QY : {:.3f}'.format(PL['relative_QY']))

    #     if uv['absorption']/uv['reference'] > 0.003:
    #         uv['absorbance_maintenance'] = min(1, math.log10(PL_uv_results['uv_time_trace'][1, -1]/uv['reference']) /math.log10(PL_uv_results['uv_time_trace'][1, 0]/uv['reference']))
    #         uv['degradation_rate'] = - math.log(uv['absorbance_maintenance'])/PL_uv_results['duration']
    #         uv['absorbance_maintenance(at_1min)'] = uv['absorbance_maintenance'] * math.exp(-(60-PL_uv_results['duration'])* uv['degradation_rate'])

    #         print('uv_absorbance_maintenance : {:.3f}%'.format(uv['absorbance_maintenance']*100))
    #         print('uv_degradation_rate: {:.3f}'.format(uv['degradation_rate']))
    #         print('uv_absorbance_maintenance(at_1min) : {:.3f}%'.format(uv['absorbance_maintenance(at_1min)']*100))
    #     else:
    #         uv['absorbance_maintenance'] = None
    #         uv['degradation_rate'] = None
    #         uv['absorbance_maintenance(at_1min)'] = None

    #     if PL['peak_time_trace'][1, -1] > 0 and PL['peak_time_trace'][1, 0] > 0:
    #         PL['maintenance'] = min(1, PL['peak_time_trace'][1, -1]/PL['peak_time_trace'][1, 0])
    #         PL['degradation_rate'] = - math.log(PL['maintenance'])/PL_uv_results['duration']
    #         PL['maintenance(at_1min)'] = PL['maintenance'] * math.exp(-(60-PL_uv_results['duration'])* PL['degradation_rate'])

    #         print('PL_maintenance : {:.3f}%'.format(PL['maintenance']*100))
    #         print('PL_degradation_rate: {:.3f}'.format(PL['degradation_rate']))
    #         print('PL_maintenance(at_1min) : {:.3f}%'.format(PL['maintenance(at_1min)']*100))
    #     else:
    #         PL['maintenance'] = None
    #         PL['degradation_rate'] = None
    #         PL['maintenance(at_1min)'] = None
    #     return uv, PL


    def plot_time_trace(self, uv_data, PL_data,  show_plot = True, save_filename = None):

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, xlabel = 'time/s', ylabel = 'Absorbance(375nm)/PL intensity')
        
        ax1.plot(uv_data['absorbance_time_trace'][0], uv_data['absorbance_time_trace'][1]/uv_data['absorbance_time_trace'][1, 0])
        ax1.plot(PL_data['peak_time_trace'][0], PL_data['peak_time_trace'][1]/PL_data['peak_time_trace'][1, 0])

        if save_filename:
            plt.savefig('%s' %save_filename)
        if show_plot:
            plt.show()
        else:
            plt.close()


    def result_plot(self, data, ylabel, text = None, xrange = [300,800], normalize = False, show_plot = True, save_filename = None):

        l_index, u_index = self._to_index_range(data[0], xrange)
        max_value = np.max(data[1][l_index:u_index+1])

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = ylabel)
        if normalize:
            data[1] = data[1]/max_value
        ax1.plot(data[0], data[1])
        plt.xlim(xrange[0],xrange[1])
        if max_value:
            if max_value > 0:
                plt.ylim(-1*max_value*0.2, max_value*1.2)
        if text:
            plt.text(xrange[0]+ (xrange[1]-xrange[0]) * 0.45, np.max(data[1])*0.9, text)
        if save_filename:
            plt.savefig('%s' %save_filename)
        if show_plot:
            plt.show()
        else:
            plt.close()


    def Abs_PL_plot(self, abs_data, PL_data, text = None, xrange = [300,800], ylabel = 'Nomalized Abs. and PL-int.', show_plot = True, save_filename = None):

        l_index, u_index = self._to_index_range(abs_data[0], xrange)

        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = ylabel)
        ax1.plot(abs_data[0], abs_data[1]/np.max(abs_data[1][l_index:u_index+1]), label = 'absorbance')
        #ax2 = ax1.twinx()
        ax1.plot(PL_data[0], PL_data[1]/np.max(PL_data[1][l_index:u_index+1]),  label = 'PL')
        plt.xlim(xrange[0],xrange[1])
        plt.ylim(-0.2, 1.2)
        plt.legend(bbox_to_anchor=(1, 1), loc='upper right')
        if text:
            plt.text(xrange[0]+ (xrange[1]-xrange[0]) * 0.6, 0.8, text)
        if save_filename:
            plt.savefig('%s' %save_filename)
        if show_plot:
            plt.show()
        else:
            plt.close()



    def _to_index_range(self, wavelength, wavelength_range):
        l_index = np.where(wavelength >= wavelength_range[0])[0][0]
        u_index = np.where(wavelength >= wavelength_range[1])[0][0]
        return l_index, u_index


    def find_max(self, data, analysis_range = [300,800]):
        l_index, u_index = self._to_index_range(data[0], analysis_range)

        d_max = np.max(data[1][l_index:u_index+1]) 
        lambda_max = data[0][np.argmax(data[1][l_index:u_index+1])+l_index] 
        return d_max, lambda_max


    def find_abs_end(self, data, threshold = 0.01, analysis_range = [300,800]):
        """
        end_index : int
            index of absortpion end. Index is in spectral range of original data
        """
        l_index, u_index = self._to_index_range(data[0], analysis_range)

        spectrum = copy.deepcopy(data)[:, l_index: u_index+1]
        spectrum[1] = spectrum[1]/np.max(spectrum[1])

        max_index = np.argmax(spectrum[1])
        end_index = None

        for i in range(len(spectrum[1]) - max_index):
            if spectrum[1][max_index + i] < threshold:
                end_index = max_index + i + l_index
                break
        if end_index:
            end_wavelength = data[0][end_index]
            return end_index, end_wavelength
        else:
            print('absorption spectrum does not go below threshold. Check the baseline')
            return None, None


    def calc_freq_spectrum(self, data, calc_range = [300,800]): 
        """ 
        calculate the PL spectrum normalized in frequency domain.
        The provided spectrum should be wavelength vs fluence (photons/s/nm)
        The dimension of the calculated spectrum :  second 
        """
        l_index, u_index = self._to_index_range(data[0], calc_range)
        spectrum = copy.deepcopy(data)[:, l_index :u_index+1]
        
        area = 0 
        for i, value in enumerate(spectrum[1]):
            if i == 0:
                area += value * (spectrum[0][i+1]-spectrum[0][i])
            elif i == len(spectrum[1])-1:
                area += value * (spectrum[0][i]-spectrum[0][i-1])
            else : 
                area += value * (spectrum[0][i+1] -spectrum[0][i-1])/2
        

        freq_spectrum = 1E9*(data[0]*1E-9)**2* data[1]/(3E8 * area)

        return np.asarray([data[0], freq_spectrum])


    def calc_gain_spectrum(self, data, calc_range = [300,800]): 
        """ 
        calculate the normalized PL spectrum in frequency domain.
        The spectrum should be wavelength vs fluence (photons/s/nm)
        The dimension of the calculated spectrum :  cm2 s
        The gain cross section can be obtained by (calculated_spectrum) * ksr/n2 
        """
        l_index, u_index = self._to_index_range(data[0], calc_range)
        spectrum = copy.deepcopy(data)[:, l_index :u_index+1]

        area = 0 
        for i, value in enumerate(spectrum[1]):
            if i == 0:
                area += value * (spectrum[0][i+1]-spectrum[0][i])
            elif i == len(spectrum[1])-1:
                area += value * (spectrum[0][i]-spectrum[0][i-1])
            else : 
                area += value * (spectrum[0][i+1] -spectrum[0][i-1])/2
        

        gain_spectrum = 1E4*1E9*(data[0]*1E-9)**4*data[1]/(8*math.pi*3E8 * area)

        return np.asarray([data[0], gain_spectrum])

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


if __name__ == '__main__':

    pl =Abs_PL(led_power=90)

    average =1000
    exposure = 0.005
    filter_size = 20

    uv_average = 30

    # pl.lamp_shutter_open()

    # time.sleep(1)

    # pl.lamp_off()

    # time.sleep(1)

    # pl.lamp_shutter_close()

    # time.sleep(1)

    # pl.lamp_on()

    # time.sleep(1)
    #pl.led_shutter_open()
    #time.sleep(3)
    #pl.led_shutter_close()

    #print(pl.led_wl)

    #print(pl.measure_uv_power(10))

    #print(pl.led_off())

    # input("please load the solvent")
    # dark = pl.measure_dark_spectrum(average, exposure, filter_size=filter_size)
    # bkg = pl.measure_transmission_spectrum(average, exposure, filter_size=filter_size, dark_correction=True, do_plot=True)
    

    # # input("please load the sample")
    # _, absorption, transimission = pl.measure_absorption_spectrum(average, exposure, bkg, filter_size=filter_size, dark_correction=True, do_plot=True)

    pl.measure_dark_spectrum(average, exposure, filter_size=filter_size, do_plot=False)
    bkg1 = pl.measure_transmission_spectrum(average, exposure, filter_size=filter_size, dark_correction=True, do_plot=False)
    time.sleep(20)
    pl.measure_dark_spectrum(average, exposure, filter_size=filter_size, do_plot=False)
    bkg2 = pl.measure_transmission_spectrum(average, exposure, filter_size=filter_size, dark_correction=True, do_plot=False)
    #time.sleep(5)
    #bkg3 = pl.measure_transmission_spectrum(average, exposure, filter_size=filter_size, dark_correction=False, do_plot=False)
    #plt.plot(bkg1[0], bkg1[1])
    #plt.plot(bkg1[0], bkg2[1])
    #plt.plot(bkg1[0], bkg3[1])
    #plt.xlim(200,800)
    #plt.show()

    plt.plot(bkg1[0], bkg2[1]/bkg1[1])
    #plt.plot(bkg1[0], bkg3[1]/bkg1[1])
    #plt.plot(bkg1[0], bkg3[1]/bkg2[1])
    plt.xlim(200,800)
    plt.show()
    #_, absorption, transimission = pl.measure_absorption_spectrum(average, exposure, bkg, filter_size=filter_size, dark_correction=False, do_plot=True)



    #pl.lamp_off()
    # print(pl.manager.spectrometer.check_status())
    # spec = pl.measure_spectrum(100, 0.01)
    # print(len(spec))


#serial_number = APT.list_available_devices()[0][1]
#print(APT.hardware_info(serial_number))

#print(serial_number)
#print(serial_number)
#if serial_number:
#Rotor = ThorlabsK10CR1(55119384, sethome=True)
#Rotor.move_by(15)
#print(Rotor.get_position())
# Rotor._move_home()
# print(Rotor._get_position())
#Rotor.move_to(10)
#print(Rotor._get_position())
