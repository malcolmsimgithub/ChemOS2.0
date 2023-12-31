# -*- coding: utf-8 -*-
"""
Zurich Instruments LabOne Python API Example.

Demonstrate how to connect to a Zurich Instruments device and record scope data
using the Scope Module when the scope is running segmented recording mode.
"""

# Copyright 2018 Zurich Instruments AG

from __future__ import print_function
import time
import numpy as np
#import zhinst.utils
from  pylab.instruments.zhinst import utils


def run_example(device_id, do_plot=False, scope_length=8192, module_historylength=1):
    """
    Run the example: Connect to a Zurich Instruments Device and obtain scope
    data using the Scope Module with the scope running in segmented recording
    mode.

    The example uses the following sequence to record multiple scope records
    (each consisting of multiple segments) using the scope module:
    0. Connect to the Data Server and configure the device.
    1. Initialize the scopeModule.
    2. Configure the scopeModule.
    3. Call the scopeModule's subscribe on the device's scope data.
    4. Loop recording multiple records:
       a. Change a measurement/device config.
       b. Call execute in the scopeModule (also resets progress to 0).
       c. Set device node scopes/0/enable to 1.
       d. Wait until a single scope record has been obtained.
       e. Set device node scopes/0/enable to 0.
       f. Call finish in the scopeModule.
       g. Call read in the scopeModule.
    5. Call scopeModule clear to stop the Module's thread and clear its memory.

    Requirements:

      MF or UHF Instrument with DIG Option (HF2 does not support segmented
        recording).

      Hardware configuration: Connect signal output 1 to signal input 1 with a
        BNC cable.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev2006` or `uhf-dev2006`.

      do_plot (bool, optional): Specify whether to plot the scope data. Default
        is no plot output.

      scope_length (int, optional): The length of the scope segment(s) to record
        (/dev..../scopes/0/length).

      module_historylength (int, optional): Value to use for the
        scopeModule/historylength parameter.

    Returns:

      data (dict of lists, each entry is a numpy array): The data returned by
      multiple calls to the scopeModule's read.

    Raises:

      Exception: If the specified device is not an MF or UHF.

      Exception: If the specified device does not have the DIG Option.

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programing Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """

    apilevel_example = 6  # The API level supported by this example.
    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
    # This example can't run with HF2 Instruments or instruments without the DIG option.
    required_devtype = r'UHF|MF'  # Regular expression of supported instruments.
    required_options = ['DIG']
    required_err_msg = "This example requires the DIG Option on either UHF or MF instruments (HF2 is unsupported)."
    (daq, device, props) = zhinst.utils.create_api_session(device_id, apilevel_example,
                                                           required_devtype=required_devtype,
                                                           required_options=required_options,
                                                           required_err_msg=required_err_msg)
    zhinst.utils.api_server_version_check(daq)

    # Enable the API's log.
    daq.setDebugLevel(3)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # Now configure the instrument for this experiment. The following channels
    # and indices work on all device configurations. The values below may be
    # changed if the instrument has multiple input/output channels and/or either
    # the Multifrequency or Multidemodulator options installed.
    # Signal output mixer amplitude [V].
    amplitude = 0.500
    out_channel = 0
    # Get the value of the instrument's default Signal Output mixer channel.
    out_mixer_channel = zhinst.utils.default_output_mixer_channel(props)
    in_channel = 0
    osc_index = 0
    scope_in_channel = 0  # scope input channel
    if props['devicetype'].startswith('UHF'):
        frequency = 1.0e6
    else:
        frequency = 100e3
    exp_setting = [
        # The output signal.
        ['/%s/sigouts/%d/on'             % (device, out_channel), 1],
        ['/%s/sigouts/%d/enables/%d'     % (device, out_channel, out_mixer_channel), 1],
        ['/%s/sigouts/%d/range'          % (device, out_channel), 1],
        ['/%s/sigouts/%d/amplitudes/%d'  % (device, out_channel, out_mixer_channel), amplitude],
        ['/%s/sigins/%d/imp50'           % (device, in_channel), 1],
        ['/%s/sigins/%d/ac'              % (device, in_channel), 0],
        ['/%s/sigins/%d/range'           % (device, in_channel), 2*amplitude],
        ['/%s/oscs/%d/freq'              % (device, osc_index), frequency]]
    node_branches = daq.listNodes('/{}/'.format(device), 0)
    if 'DEMODS' in node_branches:
        # NOTE we don't need any demodulator data for this example, but we need
        # to configure the frequency of the output signal on out_mixer_c.
        exp_setting.append(['/%s/demods/%d/oscselect' % (device, out_mixer_channel), osc_index])
    daq.set(exp_setting)

    # Perform a global synchronisation between the device and the data server:
    # Ensure that the signal input and output configuration has taken effect
    # before calculating the signal input autorange.
    daq.sync()

    # Perform an automatic adjustment of the signal inputs range based on the
    # measured input signal's amplitude measured over approximately 100 ms.
    # This is important to obtain the best bit resolution on the signal inputs
    # of the measured signal in the scope.
    zhinst.utils.sigin_autorange(daq, device, in_channel)

    # Configure the instrument's scope via the /dev..../scopes/n/ node tree branch.
    # 'length' : the length of each scope record
    daq.setInt('/%s/scopes/0/length' % device, scope_length)
    # 'channel' : select the scope channel(s) to enable.
    #  Bit-encoded as following:
    #   1 - enable scope channel 0
    #   2 - enable scope channel 1
    #   3 - enable both scope channels (requires DIG option)
    # NOTE we are only interested in one scope channel: scope_in_channel and leave
    # the other channel unconfigured
    daq.setInt('/%s/scopes/0/channel' % device, 1 << in_channel)
    # 'channels/0/bwlimit' : bandwidth limit the scope data. Enabling bandwidth
    # limiting avoids antialiasing effects due to subsampling when the scope
    # sample rate is less than the input channel's sample rate.
    #  Bool:
    #   0 - do not bandwidth limit
    #   1 - bandwidth limit
    daq.setInt('/%s/scopes/0/channels/%d/bwlimit' % (device, scope_in_channel), 1)
    # 'channels/0/inputselect' : the input channel for the scope:
    #   0 - signal input 1
    #   1 - signal input 2
    #   2, 3 - trigger 1, 2 (front)
    #   8-9 - auxiliary inputs 1-2
    #   The following inputs are additionally available with the DIG option:
    #   10-11 - oscillator phase from demodulator 3-7
    #   16-23 - demodulator 0-7 x value
    #   32-39 - demodulator 0-7 y value
    #   48-55 - demodulator 0-7 R value
    #   64-71 - demodulator 0-7 Phi value
    #   80-83 - pid 0-3 out value
    #   96-97 - boxcar 0-1
    #   112-113 - cartesian arithmetic unit 0-1
    #   128-129 - polar arithmetic unit 0-1
    #   144-147 - pid 0-3 shift value
    daq.setInt('/%s/scopes/0/channels/%d/inputselect' % (device, scope_in_channel), in_channel)
    # 'time' : timescale of the wave, sets the sampling rate to 1.8GHz/2**time.
    #   0 - sets the sampling rate to 1.8 GHz
    #   1 - sets the sampling rate to 900 MHz
    #   ...
    #   16 - sets the samptling rate to 27.5 kHz
    daq.setInt('/%s/scopes/0/time' % device, 0)
    # 'single' : only get a single scope record.
    #   0 - acquire continuous records
    #   1 - acquire a single record
    # Note: configured below in main loop.
    # daq.setInt('/%s/scopes/0/single' % device, 1)
    # Configure the scope's trigger to get aligned data
    # 'trigenable' : enable the scope's trigger (boolean).
    #   0 - acquire continuous records
    #   1 - only acquire a record when a trigger arrives
    daq.setInt('/%s/scopes/0/trigenable' % device, 1)

    # Specify the trigger channel, we choose the same as the scope input
    daq.setInt('/%s/scopes/0/trigchannel' % device, in_channel)

    # Trigger on rising edge?
    daq.setInt('/%s/scopes/0/trigrising' % device, 1)

    # Trigger on falling edge?
    daq.setInt('/%s/scopes/0/trigfalling' % device, 0)

    # Set the trigger threshold level.
    daq.setDouble('/%s/scopes/0/triglevel' % device, 0.00)

    # Set hysteresis triggering threshold to avoid triggering on noise
    # 'trighysteresis/mode' :
    #  0 - absolute, use an absolute value ('scopes/0/trighysteresis/absolute')
    #  1 - relative, use a relative value ('scopes/0trighysteresis/relative') of the trigchannel's input range
    #      (0.1=10%).
    daq.setDouble('/%s/scopes/0/trighysteresis/mode' % device, 1)
    daq.setDouble('/%s/scopes/0/trighysteresis/relative' % device, 0.05)

    # Set the trigger hold-off mode of the scope. After recording a trigger event, this specifies when the scope should
    # become re-armed and ready to trigger, 'trigholdoffmode':
    #  0 - specify a hold-off time between triggers in seconds ('scopes/0/trigholdoff'),
    #  1 - specify a number of trigger events before re-arming the scope ready to trigger ('scopes/0/trigholdcount').
    daq.setInt('/%s/scopes/0/trigholdoffmode' % device, 0)
    daq.setDouble('/%s/scopes/0/trigholdoff' % device, 50e-6)

    # Set trigdelay to 0.: Start recording from when the trigger is activated.
    daq.setDouble('/%s/scopes/0/trigdelay' % device, 0.0)

    # The trigger reference position relative within the wave, a value of 0.5 corresponds to the center of the wave.
    daq.setDouble('/%s/scopes/0/trigreference' % device, 0.25)

    # Disable trigger gating.
    daq.setInt('/%s/scopes/0/triggate/enable' % device, 0)

    # Enable segmented data transfer from the device.
    daq.setInt("/%s/scopes/0/segments/enable" % device, 1)
    # The number of segments to transfer in one shot.
    # NOTE: We will set 'segments/count' on a per-record basis below.
    # daq.setInt("/%s/scopes/0/segments/count" % device, 10)

    # Perform a global synchronisation between the device and the data server: Ensure that the settings have taken
    # effect on the device before continuing. This also clears the API's data buffers to remove any old data.
    daq.sync()

    # Check the scope_length parameter that was set:
    scope_length_set = daq.getInt('/%s/scopes/0/length' % device)
    print("Actual scope length set on the device: {} (requested {})".format(scope_length_set, scope_length))

    # Initialize and configure the Scope Module.
    scopeModule = daq.scopeModule()
    # 'scopeModule/mode' : Scope data processing mode.
    # 0 - Pass through scope segments assembled, returned unprocessed, non-interleaved.
    # 1 - Moving average, scope recording assembled, scaling applied, averaged, if averaging is enabled.
    # 2 - Not yet supported.
    # 3 - As for mode 1, except an FFT is applied to every segment of the scope recording.
    scopeModule.set('scopeModule/mode', 1)
    # 'scopeModule/averager/weight' : Average the scope shots using an exponentially weighted moving average of the
    # previous 'weight' shots.
    scopeModule.set('scopeModule/averager/weight', 1)
    # 'scopeModule/historylength' : The number of scope records to keep in the Scope Module's memory, when more records
    #   arrive in the Module from the device the oldest records are overwritten.
    scopeModule.set('scopeModule/historylength', module_historylength)

    # Subscribe to the scope's data in the module.
    wave_nodepath = '/{}/scopes/0/wave'.format(device)
    scopeModule.subscribe(wave_nodepath)

    # Loop over the desired number of measurements. For each measurement we will get a scope record consisting of
    # of the specified number of segments.
    #
    data = {}
    data[wave_nodepath] = []
    num_measurements = 5
    segment_counts = [1, 5, 10, 15, 20]
    for index, amplitude in enumerate(np.linspace(0.2, 1.0, num_measurements)):

        # Use different signal output amplitudes simply to distinguish between
        # different segments in the plot.
        daq.setDouble('/%s/sigouts/%d/amplitudes/%d'  % (device, out_channel, out_mixer_channel), amplitude)
        daq.sync()

        # Perform an automatic adjustment of the signal inputs range based on
        # the measured input signal's amplitude measured over approximately 100
        # ms.  This is important to obtain the best bit resolution on the signal
        # inputs of the measured signal in the scope.
        zhinst.utils.sigin_autorange(daq, device, in_channel)

        # Note: We should disable the scope whilst modifying settings.
        daq.setInt('/{}/scopes/0/enable'.format(device), 0)
        # Set the desired number of segments.
        daq.setInt("/{}/scopes/0/segments/count".format(device), segment_counts[index])
        daq.sync()  # Ensure the setting has taken effect on the device before continuing.
        segment_count_set = daq.getInt('/{}/scopes/0/segments/count'.format(device))
        print("Segment count set on the device: {} (requested {}).".format(segment_count_set, segment_counts[index]))

        if module_historylength == 1:
            # Set the scope to operate in 'single' mode: Once one scope record consisting of the specified number of
            # segments (>= 1) has been recorded the scope will automatically stop. Note: The device node scopes/0/single
            # will be set back to 0 by the device after recording one record.
            daq.setInt('/%s/scopes/0/single' % device, 1)
        scopeModule.set('scopeModule/clearhistory', 1)

        d = get_scope_records(device, daq, scopeModule, module_historylength)
        # Check the dictionary returned by read contains the expected data. The data returned is a dictionary with keys
        # corresponding to the recorded data's path in the node hierarchy.
        if wave_nodepath not in data:
            print("[error]: The subscribed data `{}` for measurement {} ({}) was not returned.".format(
                wave_nodepath, index, amplitude))
            continue
        else:
            num_records = len(d[wave_nodepath])
            dt = d[wave_nodepath][0][0]['dt']
            totalsamples = d[wave_nodepath][0][0]['totalsamples']
            segment_duration = dt*totalsamples/segment_counts[index]
            print("Scope data contains {} record(s).".format(num_records))
            print("Duration of each segment: {} s.".format(segment_duration))
            check_scope_record_flags(d[wave_nodepath])
            data[wave_nodepath].append(d[wave_nodepath])
        print("")

    # Stop the module's thread and remove it and its data from memory.
    scopeModule.clear()

    if do_plot and data[wave_nodepath]:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        plt.clf()
        plt.grid(True)

        clockbase = daq.getInt('/%s/clockbase' % device)

        total_segments = sum(segment_counts)
        colors = cm.rainbow(np.linspace(0, 1, total_segments))
        segment_index = 0
        for index, records in enumerate(data[wave_nodepath]):
            # We only plot the first record for each measurement. To plot all records for each measurement additionally
            # loop over `records'.
            wave = records[0][0]['wave'][scope_in_channel]
            # Reshape the array to recover the individual segments (this is only necessary in segmented mode).
            segments = wave.reshape(segment_counts[index], scope_length)
            # Create a time array relative to the trigger time.
            dt = records[0][0]['dt']
            # The timestamp is the timestamp of the last sample in the scope segment.
            timestamp = records[0][0]['timestamp']
            triggertimestamp = records[0][0]['triggertimestamp']
            t_segment = np.arange(-scope_length, 0)*dt + (timestamp - triggertimestamp)/float(clockbase)
            for segment in segments:
                plt.plot(1e3*t_segment, segment, color=colors[segment_index])
                segment_index += 1
            plt.draw()
            plt.title('{} Scope Records (consisting of different segment counts)'.format(num_measurements))
            plt.ylabel('Amplitude [V]')
            plt.xlabel('Time, relative to trigger [ms]')
        plt.axvline(0.0, linewidth=2, linestyle='--', color='k', label="Trigger time")
        plt.autoscale(enable=True, axis='x', tight=True)

    return data


def get_scope_records(device, daq, scopeModule, num_records=1):
    """
    Obtain scope records from the device using an instance of the Scope Module.
    """

    # Tell the module to be ready to acquire data; reset the module's progress to 0.0.
    scopeModule.execute()

    # Enable the scope: Now the scope is ready to record data upon receiving triggers.
    daq.setInt('/%s/scopes/0/enable' % device, 1)
    daq.sync()

    start = time.time()
    timeout = 30  # [s]
    records = 0
    progress = 0
    # Wait until the Scope Module has received and processed the desired number of records.
    while (records < num_records) or (progress < 1.0):
        time.sleep(0.5)
        records = scopeModule.getInt("scopeModule/records")
        progress = scopeModule.progress()[0]
        print(("Scope module has acquired {} records (requested {}). "
               "Progress of current segment {}%.").format(records, num_records, 100.*progress), end='\r')
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
            print("\nScope Module did not return {} records after {} s - forcing stop.".format(num_records, timeout))
            break
    print("")
    daq.setInt('/%s/scopes/0/enable' % device, 0)

    # Read out the scope data from the module.
    data = scopeModule.read(True)

    # Stop the module; to use it again we need to call execute().
    scopeModule.finish()

    return data


def check_scope_record_flags(scope_records):
    """
    Loop over all records and print a warning to the console if an error bit in
    flags has been set.

    Warning: This function is intended as a helper function for the API's
    examples and it's signature or implementation may change in future releases.
    """
    num_records = len(scope_records)
    for index, record in enumerate(scope_records):
        if record[0]['flags'] & 1:
            print('Warning: Scope record {}/{} flag indicates dataloss.'.format(index, num_records))
        if record[0]['flags'] & 2:
            print('Warning: Scope record {}/{} indicates missed trigger.'.format(index, num_records))
        if record[0]['flags'] & 4:
            print('Warning: Scope record {}/{} indicates transfer failure (corrupt data).'.format(index, num_records))
        totalsamples = record[0]['totalsamples']
        for wave in record[0]['wave']:
            # Check that the wave in each scope channel contains the expected number of samples.
            assert len(wave) == totalsamples, \
                'Scope record {}/{} size does not match totalsamples.'.format(index, num_records)
