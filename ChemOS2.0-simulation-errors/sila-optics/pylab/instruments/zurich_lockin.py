from pylab.instruments.zhinst import utils as zutils 
from pylab.instruments import zhinst
# from pylab.instruments.zhinst import ziPython
try:
    from pylab.instruments.zhinst.ziPython import ziDiscovery
except:
    pass
#from pylab.instruments.zhinst.ziPython import ziListEnum
import numpy as np
import time
import math
import matplotlib.pyplot as plt

#import zhinst.examples.common.example_scope as example
#example.run_example('dev4357', do_plot = True)

def find_device():
    """  
    get a list of available device id

    Returns
    --------------------
    id : list
        list of id of available device
    """
    zd = ziDiscovery()
    lockin_id = zd.findAll()
    print(lockin_id)


class MFLI500(object):
    """  
    A Class control to the Zurich MFLI500 lock-in amplifier.

    Parameters
    --------------------
    lockin_id : str
    """
    
    def __init__(self, lockin_id):
        
        self.api_level = 6
        self.devtype = 'MFLI'

        #create a api session
        (self.daq, self.device, self.props) = zutils.create_api_session(lockin_id, self.api_level, 
                          required_devtype = self.devtype, 
                          required_options=None, 
                          required_err_msg = 'error')
        
        print(self.device)

        # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
        zutils.disable_everything(self.daq, self.device)


    def _enum_channel(self, channel):
        channels = {'Vin' : 0, 
                    'Iin' : 1,
                    'trigin1' : 2,
                    'trigin2': 3,
                    'augout1' : 4,     
                    'augout2' : 5,    
                    'augout3' : 6,    
                    'augout4' : 7,    
                    'augin1' : 8,
                    'augin2' : 9,
                    'sigout1' : 12,
                    'trigout1' : 14,
                    'trigout2' : 15
                    }
        return channels[channel]

    
    def ask_timeconstant(self):
        """  
        get a time constant setting of the demodulator

        Returns
        --------------------
        time constant : float
        """
        time_constant = self.daq.get('/%s/demods/0/timeconstant' %self.device, True)['/%s/demods/0/timeconstant' %self.device]['value']
        #print(time_constant)
        return time_constant


    def set_param(self, setting):
        """  
        set the device paramter

        Parameters
        --------------------
        setting : list
            list of the setting 
        """
        time_constant = self.ask_timeconstant()

        self.daq.set(setting)

        # Wait for the demodulator filter to settle.
        time.sleep(10*time_constant)
 
        # Perform a global synchronisation between the device and the data server:
        # Ensure that the settings have taken effect on the device before issuing
        # the getSample() command. Note: the sync() must be issued after waiting for
        # the demodulator filter to settle above.
        self.daq.sync()


    def ask_param_list(self, path):
        """  
        get a list of parameter in a specific path

        Parameters
        --------------------
        path : str
            e.g. '/dev4357/auxouts/0/outputselect'

        Returns
        --------------------
        list of the parameter : list
        """
        return self.daq.getList(path)


    def set_auxout(self, channel, source, preoffset, scale, offset, llim = -10, ulim = 10):
        """  
        set an auxiliary output. The output voltage is calcualted as;
        V = (source + preoffset)*scale + offset
        
        Parameters
        --------------------
        channel : int
            0-4
        source : int
            -1:manual, 0:DemodX, 1:DemodY, 2:DemodR, 3:Demodq,...
        preoffset : float
        scale : float
        llim : float
            lower limit of the output voltage
        ulim 
            upper limit of the output voltage
        """

        time_constant = self.ask_timeconstant()

        auxout_setting = [['/%s/auxouts/%d/outputselect'         % (self.device, channel), source], 
                            ['/%s/auxouts/%d/preoffset'          % (self.device, channel), preoffset],
                            ['/%s/auxouts/%d/scale'              % (self.device, channel), scale],
                            ['/%s/auxouts/%d/offset'             % (self.device, channel), offset],
                            ['/%s/auxouts/%d/limitlower'         % (self.device, channel), llim],
                            ['/%s/auxouts/%d/limitupper'         % (self.device, channel), ulim],
                         ]
        
        self.daq.set(auxout_setting)

        time.sleep(10*time_constant)

        self.daq.sync()


    def set_scope(self, sigin = 'Vin', sigrange = 'auto', sampling = 8, length = 16384, sigfl = False, sigac = False, sigdiff = False, imp50 = False, trig_enable= True):
        """  
        Set scope parameters.
        
        Parameters
        --------------------
        sigin : str
            input channel ('Vin', 'Iin', 'trigin1~2', augout1~4', 'augin1~2', 'sigout1', 'trigout1~2')
        sigrange : str
            3.0 mV, 10 mV, 30 mV, 100 mV, 300 mV, 1.0 V, 3.0 V, or auto
        sampling : int
            sampling rate = 60MHz*2**(-sampling)
        length : int
            length of the data to be recorded (min 4096)
        sigfl : bool
            switching between floating and connected to ground
        sigdiff : bool
            enable/disable differential input mode
        imp50 : bool
            Switches between 50 Ω (True) and 10 MΩ (False).
        trig_enable : bool
            enable/disable trigger
        """

        in_channel = self._enum_channel(sigin)

        exp_setting = [['/%s/sigouts/%d/on'             % (self.device, 0), 0],
                        ['/%s/sigouts/%d/enables/%d'     % (self.device, 0, 1), 0],
                        ['/%s/sigins/0/imp50'           % self.device, imp50],
                        ['/%s/sigins/0/ac'              % self.device, sigac],
                        ['/%s/sigins/0/diff'              % self.device,  sigdiff],
                        ['/%s/sigins/0/float'              % self.device, sigfl]
                        ]
        if sigrange == 'auto':
            exp_setting.append(['/%s/sigins/%d/autorange'  %(self.device, in_channel), 1])
        else: 
            exp_setting.append(['/%s/sigins/%d/autorange'  %(self.device, in_channel), 0])
            exp_setting.append(['/%s/sigins/%d/range'  %(self.device, in_channel), sigrange])

        self.daq.set(exp_setting)
        self.daq.setInt('/%s/scopes/0/trigenable' %self.device, trig_enable)
        self.daq.setInt('/%s/scopes/0/time' %self.device, sampling)  
        self.daq.setInt('/%s/scopes/0/length' %self.device, length)  

        self.daq.sync()       


    def set_trigger(self, channel = 'trigin1', mode = 'rise', level = 0.5, hyst = 0.1, ref = 0.5, delay = 0.0, holdoff = 2, enable_gate = False):
        """  
        Set trigger parameters.
        
        Parameters
        --------------------
        channel : str
            channel used for the trigger 
            ('Vin', 'Iin', 'trigin1~2', augout1~4', 'augin1~2', 'sigout1', 'trigout1~2')
        mode : str
            trigger mode ('rise', 'fall', 'both')
        level : float
            trigger level (V)
        hyst : float 
            The source signal must deviate from the trigger level before the trigger is rearmed again. 
            The vale is relative to the adjusted full scale signal input range. Set to 0 to turn it off. 
        ref : float
            The trigger reference position relative within the wave, a value of 0.5 corresponds to the center of the wave.
        delay : float
            Trigger position relative to reference. A positive delay results in less data being acquired before the trigger point, 
            a negative delay results in more data being acquired before the trigger point. (s)
        holdoff : int
            number of the trigger event that will trigger the next recording after a recording event. 
            A value one will start a recording for each trigger event.
        enable_gate : bool
            If enabled the trigger will be gated by the trigger gating input signal. This feature requires the MF-DIG option.
        """

        trig_channel = self._enum_channel(channel)

        self.daq.setInt('/%s/scopes/0/trigchannel' % self.device, trig_channel)
        if mode == 'rise':
            self.daq.setInt('/%s/scopes/0/trigslope' % self.device, 1)
        elif mode == 'fall':
            self.daq.setInt('/%s/scopes/0/trigslope' % self.device, 2)
        elif mode == 'both':
            self.daq.setInt('/%s/scopes/0/trigslope' % self.device, 3)

        self.daq.setInt('/%s/scopes/0/triggate/enable' % self.device, enable_gate)

        self.daq.setDouble('/%s/scopes/0/triglevel' % self.device, level)
        self.daq.setDouble('/%s/scopes/0/trigdelay' % self.device, delay)

        # Set hysteresis triggering threshold to avoid triggering on noise
        # 'trighysteresis/mode' :
        #  0 - absolute, use an absolute value ('scopes/0/trighysteresis/absolute')
        #  1 - relative, use a relative value ('scopes/0trighysteresis/relative') of the trigchannel's input range
        #      (0.1=10%).
        self.daq.setDouble('/%s/scopes/0/trighysteresis/mode' % self.device, 1)
        self.daq.setDouble('/%s/scopes/0/trighysteresis/relative' % self.device, hyst)  # 0.1=10%

        # Set the trigger hold-off mode of the scope. After recording a trigger event, this specifies when the scope should
        # become re-armed and ready to trigger, 'trigholdoffmode':
        #  0 - specify a hold-off time between triggers in seconds ('scopes/0/trigholdoff'),
        #  1 - specify a number of trigger events before re-arming the scope ready to trigger ('scopes/0/trigholdcount').
        self.daq.setInt('/%s/scopes/0/trigholdoffmode' % self.device, 1)
        self.daq.setDouble('/%s/scopes/0/trigholdoffcount' % self.device, holdoff)

        self.daq.setDouble('/%s/scopes/0/trigreference' % self.device, ref)

        self.daq.setInt('/%s/scopes/0/trigenable' %self.device, 1)

        self.daq.sync()


    def get_scope_records(self, mode = 1, averager = 1, record_num = 1, history_length = 1, timeout = 60):
        """
        Obtain scope records from the device using an instance of the Scope Module.
         
        Parameters
        --------------------
        mode : int
            scope mode (1:time_domain, 3 : fft_domain) 
        averager : int 
            number of record averaged (1:none, 1>:exponentially weighted moving average)
        record_num : int
            number of records to be acquired
        history_length : int
            The number of scope records to keep in the Scope Module's memory. 
            When more records arrive in the Module from the device the oldest records are overwritten.
        Returns
        --------------------
        record : dict
            {time (or frequency) : list,  wave : list of list}
            recorded signals (record_num)
        """
        scopeModule = self.daq.scopeModule()
        scopeModule.set('scopeModule/mode', mode)
        scopeModule.set('scopeModule/averager/weight', averager)
        print("MFLI:averager:: %s" %scopeModule.getInt("scopeModule/averager/weight"))
        print('MFLI:record_num:: %s' %record_num)
        print('MFLI:histry_length:: %s' %history_length)
        scopeModule.set('scopeModule/historylength', history_length)

        # Subscribe to the scope's data in the module.
        wave_nodepath = '/%s/scopes/0/wave' %self.device
        scopeModule.subscribe(wave_nodepath)
        # Tell the module to be ready to acquire data; reset the module's progress to 0.0.
        scopeModule.execute()

        # Enable the scope: Now the scope is ready to record data upon receiving triggers.
        self.daq.setInt('/%s/scopes/0/enable' % self.device, 1)
        self.daq.sync()

        start = time.time()
        records = 0
        progress = 0
        # Wait until the Scope Module has received and processed the desired number of records.
        while (records < record_num) or (progress < 1.0):
            time.sleep(0.5)
            records = scopeModule.getInt("scopeModule/records")
            progress = scopeModule.progress()[0]
            print(("Scope module has acquired {} records (requested {}). "
                "Progress of current segment {}%.").format(records, record_num, 100.*progress), end='\r')
            # Advanced use: It's possible to read-out data before all records have been recorded (or even before all
            # segments in a multi-segment record have been recorded). Note that complete records are removed from the Scope
            # Module and can not be read out again; the read-out data must be managed by the client code. If a multi-segment
            # record is read-out before all segments have been recorded, the wave data has the same size as the complete
            # data and scope data points currently unacquired segments are equal to 0.
            #
            # data = scopeModule.read(True)
            # wave_nodepath = '/{}/scopes/0/wave'.format(device)
            # if wave_nodepath in data:
            #   Do something with the data...
            if (time.time() - start) > timeout:
                # Break out of the loop if for some reason we're no longer receiving scope data from the device.
                print("\nScope Module did not return {} records after {} s - forcing stop.".format(record_num, timeout))
                break
        print("")
        self.daq.setInt('/%s/scopes/0/enable' % self.device, 0)

        # Read out the scope data from the module.
        data = scopeModule.read(True)

        # Stop the module; to use it again we need to call execute().
        scopeModule.finish()
        
        result = {'wave': []}

        clockbase = 6e7/2**(self.daq.getInt('/%s/scopes/0/time' %self.device))
        #print(len(data[wave_nodepath]), len(data[wave_nodepath][0]))

        for record in data[wave_nodepath]:
#        record = data[wave_nodepath][:][0]
            wave = record[0]['wave'][0,:]
#        print(len(record['wave']))
            result['wave'].append(np.asarray(wave))

        totalsamples = data[wave_nodepath][-1][0]['totalsamples']

        if mode == 1:
            dt = data[wave_nodepath][-1][0]['dt']
            timestamp = data[wave_nodepath][-1][0]['timestamp']
            triggertimestamp = data[wave_nodepath][-1][0]['triggertimestamp']
            t = np.arange(-totalsamples, 0)*dt + (timestamp - triggertimestamp)/float(clockbase)
            result['time'] = t
        elif mode == 3: 
            # We're in FFT mode.
            scope_time = 0
            scope_rate = clockbase/2**scope_time
            f = np.linspace(0, scope_rate/2, totalsamples)
            result['frequency'] = f
        return result



    def set_lockin(self, sigin = 'Vin', sigAC = False, sigdiff = False, sigfl = False, sigrange = 'auto', sigscale = 1,
                         ref = 'trigin1', ref_harm = 1,  LP_order = 4, LP_BW = 1,  sinc = False, enable = True):
        """  
        Set lockin parameters.
        
        Parameters
        --------------------
        sigin : str
            input channel ('Vin', 'Iin', 'trigin1~2', augout1~4', 'augin1~2', 'sigout1', 'trigout1~2')
        sigAC: bool
            enable/disable AC mode (can be used to remove DC offset)
        sigdiff : bool
            enable/disable differential input mode   
        sigfl : bool
            switching between floating and connected to ground 
        sigrange : str
            3.0 mV, 10 mV, 30 mV, 100 mV, 300 mV, 1 V, 3.0 V, or auto
        sigscale : float
            Applies an arbitrary scale factor to the input signal.
        ref : str or int
            external reference source(str) or internal reference frequency(int in Hz) 
        ref_harm : int
            Divides the demodulator's reference frequency by an integer factor in external reference mode. 
        LP_order : int
            order of the low-pass filter
        LP_BW : float
            bandwidth of the low-path filter (Hz)
        sinc : bool
            enable/disable sinc filter
        enable : bool
            enable/disable the data acquisition
        """
        time_constant = 1/((LP_order+1)*math.pi*LP_BW)

        #To do calculation of timeconstant
        exp_setting = [['/%s/demods/0/adcselect'        % self.device, self._enum_channel(sigin)],
                         ['/%s/sigins/0/ac'             % self.device, sigAC],
                         ['/%s/sigins/0/diff'           % self.device, sigdiff],
                         ['/%s/sigins/0/float'          % self.device, sigfl],
                         ['/%s/sigins/0/scaling'        % self.device, sigscale],
                         ['/%s/demods/0/enable'         % self.device, 1],
                         ['/%s/demods/0/rate'           % self.device, LP_BW*10],
                         ['/%s/demods/0/order'          % self.device, LP_order],
                         ['/%s/demods/0/timeconstant'   % self.device, time_constant],
                         ['/%s/demods/0/oscselect'      % self.device, 0],
                         ['/%s/demods/0/harmonic'       % self.device, ref_harm],
                         ['/%s/demods/0/sinc'           % self.device, sinc],
                         ['/%s/sigouts/0/on'            % self.device, 0]
                         ]

        if sigrange == 'auto':
            exp_setting.append(['/%s/sigins/0/autorange'  %self.device, 1])
        else: 
            exp_setting.append(['/%s/sigins/0/autorange'  %self.device, 0])
            exp_setting.append(['/%s/sigins/0/range'  %self.device, sigrange])

        #To do : select external source
        if type(ref) == str:
            exp_setting.append(['/%s/extrefs/0/enable'  % self.device, 1])
            exp_setting.append(['/%s/demods/1/adcselect' %self.device, self._enum_channel(ref)])
        else: 
            exp_setting.append(['/%s/extrefs/0/enable'  % self.device, 0])
            exp_setting.append(['/%s/oscs/0/freq'    % self.device, ref])
        
        if enable == True:
            exp_setting.append(['/%s/demods/0/enable'  %self.device, 1])
        else:
            exp_setting.append(['/%s/demods/0/enable'  %self.device, 0])

        self.daq.set(exp_setting)
        
        # Unsubscribe any streaming data.
        self.daq.unsubscribe('*')
        data = {}
    
        # Wait for the demodulator filter to settle.
        settle_time = 50 * self.ask_timeconstant()
        time.sleep(settle_time)
        # Perform a global synchronisation between the device and the data server:
        # Ensure that 1. the settings have taken effect on the device before issuing
        # the poll() command and 2. clear the API's data buffers. Note: the sync()
        # must be issued after waiting for the demodulator filter to settle above.
        self.daq.sync()



    def measure(self, sampling_rate, duration, burst_duration = 1, do_plot=False, filename = None):
        """
        Obtain lockin records.
         
        Parameters
        --------------------
        sampling_rate : int
            sampling rate in Hz, 
        duration : float
            duration for the data aquisition in seconds 
        burst_duration : float
            Time in seconds for each data burst/segment.
        do_plot : bool
            plot result if True
        filename : str or None

        Returns
        --------------------
        data : 2D list
        """
        #sampling_rate in Hz, duration in s
        # The list of signal paths that we would like to record in the module.
        demod_path = '/{}/demods/0/sample'.format(self.device)
        self.daq.set([['/{}/demods/0/enable'.format(self.device), 1]]) #enable aquisition

        signal_paths = []
        signal_paths.append(demod_path + '.r')  # The demodulator X output.
        signal_paths.append(demod_path + '.theta')  # The demodulator Y output.
        
        # Check the device has demodulators.
        flags = zhinst.ziPython.ziListEnum.recursive | zhinst.ziPython.ziListEnum.absolute | zhinst.ziPython.ziListEnum.streamingonly
        streaming_nodes = self.daq.listNodes('/{}'.format(self.device), flags)
        if demod_path.upper() not in streaming_nodes:
            print("Device {} does not have demodulators. Please modify the example to specify".format(self.device),
                  "a valid signal_path based on one or more of the following streaming nodes: ",
                  "{}".format('\n'.join(streaming_nodes)))
            raise Exception("Demodulator streaming nodes unavailable - see the message above for more information.")
    
        total_duration = duration   #time in second
        module_sampling_rate = sampling_rate  # Number of points/second
        num_cols = int(np.ceil(module_sampling_rate*burst_duration))
        num_bursts = int(np.ceil(total_duration/burst_duration))
    
        # Create an instance of the Data Acquisition Module.
        h = self.daq.dataAcquisitionModule()
        h.set("dataAcquisitionModule/device", self.device)
        h.set("dataAcquisitionModule/type", 0) #Trigger off(continuous aquisition)
        h.set("dataAcquisitionModule/grid/mode", 2) #Linear interpolation mode
        h.set("dataAcquisitionModule/count", num_bursts) #the number of bursts of data the module should return
        h.set("dataAcquisitionModule/duration", burst_duration)  #Burst duration in seconds.
        h.set("dataAcquisitionModule/grid/cols", num_cols) #The number of points within each duration.
    
        if filename is not None:
            h.set('dataAcquisitionModule/save/fileformat', 1) # 0:MatLab, 1:CSV
            h.set('dataAcquisitionModule/save/filename', filename)
            h.set('dataAcquisitionModule/save/saveonread', 1) #enables autosave
    
        data = [[],[],[]] 
        for signal_path in signal_paths:
            print("Subscribing to", signal_path)
            h.subscribe(signal_path)  
    
        clockbase = float(self.daq.getInt("/{}/clockbase".format(self.device)))
    
        ts0 = np.nan
        read_count = 0
    
        def read_data_update_plot(data, timestamp0):
            """
            Read the acquired data out from the module and plot it. Raise an
            AssertionError if no data is returned.
            """
            data_read = h.read(True)
            returned_signal_paths = [signal_path.lower() for signal_path in data_read.keys()]
            progress = h.progress()[0]
            # Loop over all the subscribed signals:
            for i, signal_path in enumerate(signal_paths):
                if signal_path.lower() in returned_signal_paths:
                    for index, signal_burst in enumerate(data_read[signal_path.lower()]):
                        if np.any(np.isnan(timestamp0)):
                            # Set our first timestamp to the first timestamp we obtain
                            for j, signal in enumerate(signal_burst['value'][0, :]):
                                if np.isnan(signal) == False :
                                    timestamp0 = signal_burst['timestamp'][0, j]
                                    break
                        # Convert from device ticks to time in seconds.
                        t = (signal_burst['timestamp'][0, :] - timestamp0)/clockbase
                        value = signal_burst['value'][0, :]
                        num_samples = len(signal_burst['value'][0, :])
                        dt = (signal_burst['timestamp'][0, -1] - signal_burst['timestamp'][0, 0])/clockbase
                        if i ==0:
                            data[0].extend(t)
                        data[i + 1].extend(value)
                else:
                    pass
            
            return data, timestamp0
    
        # Start recording data.
        h.execute()
    
        # Record data in a loop with timeout.
        timeout = 3*total_duration
        t0_measurement = time.time()
        # The maximum time to wait before reading out new data.
        t_update = 0.9*burst_duration
        while not h.finished():
            t0_loop = time.time()
            if time.time() - t0_measurement > timeout:
                raise Exception("Timeout after {} s - recording not complete. Are the streaming nodes enabled? "
                                "Has a valid signal_path been specified?".format(timeout))
            data, ts0 = read_data_update_plot(data, ts0)
            read_count += 1
            # We don't need to update too quickly.
            time.sleep(max(0, t_update - (time.time() - t0_loop)))
    
        # There may be new data between the last read() and calling finished().
        data, _ = read_data_update_plot(data, ts0)
    
        # Before exiting, make sure that saving to file is complete (it's done in the background)
        # by testing the 'dataAcquisitionModule/save/save' parameter.
        timeout = 1.5*total_duration
        t0 = time.time()
        while h.getInt('dataAcquisitionModule/save/save') != 0:
            time.sleep(0.1)
            if time.time() - t0 > timeout:
                raise Exception("Timeout after {} s before data save completed.".format(timeout))
    
        h.clear() # destroy the instance in module

        if do_plot:
            plt.plot(data[0], data[1])
            plt.xlabel('time/s')
            plt.ylabel('demodulator R/mV')
            plt.show()
        return data


    def poll(self):
        """
        Obtain demodulator data using ziDAQServer's blocking (synchronous) poll() command
        
        Returns
        --------------------
        sample : dict
            {time : list, R : list, phi : list}
        """
       # Unsubscribe any streaming data.
        self.daq.unsubscribe('*')
    
        # Subscribe to the demodulator's sample node path.
        path = '/%s/demods/0/sample' % self.device
        self.daq.subscribe(path)
    
        # Sleep for demonstration purposes: Allow data to accumulate in the data
        # server's buffers for one second: poll() will not only return the data
        # accumulated during the specified poll_length, but also for data
        # accumulated since the subscribe() or the previous poll.
        # For demonstration only: We could, for example, be processing the data
        # returned from a previous poll().
        sleep_length = 3.0

        time.sleep(sleep_length)
    
        # Poll the subscribed data from the data server. Poll will block and record
        # for poll_length seconds.
        poll_length = 0.1  # [s]
        poll_timeout = 500  # [ms]
        poll_flags = 0
        poll_return_flat_dict = True
        data = self.daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)
    
        # Unsubscribe from all paths.
        self.daq.unsubscribe('*')
    
        # Check the dictionary returned is  non-empty
        assert data, "poll() returned an empty data dictionary, did you subscribe to any paths?"
    
        # The data returned is a dictionary of dictionaries that reflects the node's path.
        # Note, the data could be empty if no data had arrived, e.g., if the demods
        # were disabled or had demodulator rate 0.
        assert path in data, "The data dictionary returned by poll has no key `%s`." % path
    
        # Access the demodulator sample using the node's path.
        sample = data[path]
    
        # Let's check how many seconds of demodulator data were returned by poll.
        # First, get the sampling rate of the device's ADCs, the device clockbase...
        clockbase = float(self.daq.getInt('/%s/clockbase' % self.device))
        # ... and use it to convert sample timestamp ticks to seconds:
        dt_seconds = (sample['timestamp'][-1] - sample['timestamp'][0])/clockbase
        print("poll() returned {:.3f} seconds of demodulator data.".format(dt_seconds))
        tol_percent = 10
        dt_seconds_expected = sleep_length + poll_length
        assert (dt_seconds - dt_seconds_expected)/dt_seconds_expected*100 < tol_percent, \
            "Duration of demod data returned by poll() (%.3f s) differs " % dt_seconds + \
            "from the expected duration (%.3f s) by more than %0.2f %%." % \
            (dt_seconds_expected, tol_percent)
    
        # Calculate the demodulator's magnitude and phase and add them to the dict.
        sample['R'] = np.abs(sample['x'] + 1j*sample['y'])
        sample['phi'] = np.angle(sample['x'] + 1j*sample['y'])
        print("Average measured RMS amplitude is {:.3e} V.".format(np.mean(sample['R'])))
    
            # Convert timestamps from ticks to seconds via clockbase.
        sample['time'] = (sample['timestamp'] - sample['timestamp'][0])/clockbase
    
        return sample

