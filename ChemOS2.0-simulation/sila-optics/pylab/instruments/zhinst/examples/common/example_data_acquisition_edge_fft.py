# -*- coding: utf-8 -*-
"""
Zurich Instruments LabOne Python API Example

Python API Example for using the FFT feature of the Data Acquisition Core Module. This
example demonstrates how to get FFTs from triggered bursts of demodulator data when a
demodulator's R value is larger than a specified threshold using an edge
trigger.

The Data Acquisition Module implements software triggering which operates
analogously to the types of triggering found in laboratory oscilloscopes. The
Data Acquisition Module has a non-blocking (asynchronous) interface, it starts it's
own thread to communicate with the data server.

Note: This example requires a feedback cable between Signal Output 1 and Signal
Input 1 and changes the signal output's amplitude in order to create a signal
upon which to trigger.

"""

# Copyright 2016 Zurich Instruments AG

from __future__ import print_function
import time
import numpy as np
from  pylab.instruments.zhinst import utils

def run_example(device_id, amplitude=0.25, do_plot=False):
    """Run the example: Record FFTs of bursts of demodulator sample data when a
    demodulator's R becomes larger than a specified threshold using ziPython's
    Data Acquisition Module.

    Requirements:

      Hardware configuration: Connect signal output 1 to signal input 1 with a
      BNC cable.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev2006` or `uhf-dev2006`.

      amplitude (float, optional): The amplitude to set on the signal output.

      do_plot (bool, optional): Specify whether to plot the recorded
        data. Default is no plot output.

    Returns:

      data (dict): The data structure as returnd by the Data Acuisition Module's
        read() command. The difrom pylab.instrumentsctionary contains the values of the module's
        parameters and the data returned by the module.

    Raises:

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
    err_msg = "This example only supports instruments with demodulators."
    (daq, device, props) = utils.create_api_session(device_id, apilevel_example,
                                                           required_devtype='.*LI|.*IA|.*IS',
                                                           required_err_msg=err_msg)
    utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    utils.disable_everything(daq, device)

    # Now configure the instrument for this experiment. The following channels
    # and indices work on all device configurations. The values below may be
    # changed if the instrument has multiple input/output channels and/or either
    # the Multifrequency or Multidemodulator options installed.
    out_channel = 0
    out_mixer_channel = utils.default_output_mixer_channel(props)
    in_channel = 0
    demod_index = 0
    osc_index = 0
    demod_rate = 10e3
    time_constant = 8e-5
    frequency = 400e3
    exp_setting = [['/%s/sigins/%d/ac'             % (device, in_channel), 0],
                   ['/%s/sigins/%d/imp50'          % (device, in_channel), 0],
                   ['/%s/sigins/%d/range'          % (device, in_channel), 3*amplitude],
                   ['/%s/demods/%d/enable'         % (device, demod_index), 1],
                   ['/%s/demods/%d/rate'           % (device, demod_index), demod_rate],
                   ['/%s/demods/%d/adcselect'      % (device, demod_index), in_channel],
                   ['/%s/demods/%d/order'          % (device, demod_index), 4],
                   ['/%s/demods/%d/timeconstant'   % (device, demod_index), time_constant],
                   ['/%s/demods/%d/oscselect'      % (device, demod_index), osc_index],
                   ['/%s/demods/%d/harmonic'       % (device, demod_index), 1],
                   ['/%s/oscs/%d/freq'             % (device, osc_index), frequency],
                   ['/%s/sigouts/%d/on'            % (device, out_channel), 1],
                   ['/%s/sigouts/%d/enables/%d'    % (device, out_channel, out_mixer_channel), 1],
                   ['/%s/sigouts/%d/range'         % (device, out_channel), 1],
                   ['/%s/sigouts/%d/amplitudes/%d' % (device, out_channel, out_mixer_channel), amplitude]]
    daq.set(exp_setting)

    # Wait for the demodulator filter to settle.
    time.sleep(10*time_constant)

    # Perform a global synchronisation between the device and the data server:
    # Ensure that 1. the settings have taken effect on the device before issuing
    # the poll() command and 2. clear the API's data buffers. Note: the sync()
    # must be issued after waiting for the demodulator filter to settle above.
    daq.sync()

    # Create an instance of the Data Acquisition Module.
    trigger = daq.dataAcquisitionModule()

    # Below we will generate num_pulses pulses on the signal outputs in order to
    # demonstrate the triggering functionality. We'll configure the Software
    # Trigger's threshold according to these levels.
    sigouts_high = 1.5*amplitude
    sigouts_low = 1.0*amplitude
    num_pulses = 20

    # Configure the Data Acquisition Module.
    # Set the device that will be used for the trigger - this parameter must be set.
    trigger.set('dataAcquisitionModule/device', device)
    # We will trigger on the demodulator sample's R value.
    trigger_path = '/%s/demods/%d/sample.r' % (device, demod_index)
    triggernode = trigger_path
    trigger.set('dataAcquisitionModule/triggernode', triggernode)
    # Use an edge trigger.
    trigger.set('dataAcquisitionModule/type', 1)  # 1 = edge
    # Trigger on the positive edge.
    trigger.set('dataAcquisitionModule/edge', 1)  # 1 = positive
    # The set the trigger level.
    # Scale by 1/sqrt(2) due to the demodulator's R RMS value.
    trigger_level = 0.5*(sigouts_low + sigouts_high)/np.sqrt(2)
    print("Setting trigger/0/level to {:.3f}.".format(trigger_level))
    trigger.set('dataAcquisitionModule/level', trigger_level)
    # Set the trigger hysteresis to a percentage of the trigger level: This
    # ensures that triggering is robust in the presence of noise. The trigger
    # becomes armed when the signal passes through the hysteresis value and will
    # then actually trigger if the signal additionally passes through the
    # trigger level. The hysteresis value is applied negatively for a positive
    # edge trigger relative to the trigger level (positively for a negative edge
    # trigger).
    trigger_hysteresis = 0.05*trigger_level
    print("Setting trigger/0/hysteresis {:.3f}.".format(trigger_hysteresis))
    trigger.set('dataAcquisitionModule/hysteresis', trigger_hysteresis)
    # The number of times to trigger.
    trigger_count = int(num_pulses/2)
    trigger.set('dataAcquisitionModule/count', trigger_count)
    trigger.set('dataAcquisitionModule/holdoff/count', 0)
    trigger.set('dataAcquisitionModule/holdoff/time', 0.100)
    trigger_delay = -0.020
    trigger.set('dataAcquisitionModule/delay', trigger_delay)
    demod_rate = daq.getDouble('/%s/demods/%d/rate' % (device, demod_index))
    # 'dataAcquisitionModule/grid/mode' - Specify the interpolation method of
    #   the returned data samples.
    #
    # 1 = Nearest. If the interval between samples on the grid does not match
    #     the interval between samples sent from the device exactly, the nearest
    #     sample (in time) is taken.
    #
    # 2 = Linear interpolation. If the interval between samples on the grid does
    #     not match the interval between samples sent from the device exactly,
    #     linear interpolation is performed between the two neighbouring
    #     samples.
    #
    # 4 = Exact. The subscribed signal with the highest sampling rate (as sent
    #     from the device) defines the interval between samples on the DAQ
    #     Module's grid. If multiple signals are subscribed, these are
    #     interpolated onto the grid (defined by the signal with the highest
    #     rate, "highest_rate"). In this mode, dataAcquisitionModule/duration is
    #     read-only and is defined as num_cols/highest_rate.
    trigger.set('dataAcquisitionModule/grid/mode', 4)
    # For an FFT the number of samples needs to be a binary power
    # sample_count = int(demod_rate * trigger_duration)
    sample_count = 2048
    # The duration (the length of time to record each time we trigger) must fit exactly with
    # the number of samples. Otherwise in exact mode, it will be adjusted to fit.
    trigger_duration = sample_count/demod_rate
    trigger.set('dataAcquisitionModule/duration', trigger_duration)
    trigger.set('dataAcquisitionModule/grid/cols', sample_count)
    trigger_duration = trigger.getDouble('dataAcquisitionModule/duration')
    # The size of the internal buffer used to record triggers (in seconds), this
    # should be larger than trigger_duration.
    buffer_size = trigger.getInt('dataAcquisitionModule/buffersize')

    # In this example we obtain the absolute values of the FFT of the quadrature components
    # (x+iy) of the recorded data.
    # This is done by appending the signal/operations to the basic node path using dot notation and
    # then subscribing to this path (signal).
    # We could additionally subscribe to other node paths.
    signal_path = '/%s/demods/%d/sample.xiy.fft.abs' % (device, demod_index)
    filter_path = '/%s/demods/%d/sample.xiy.fft.abs.filter' % (device, demod_index)
    trigger.subscribe(signal_path)
    trigger.subscribe(filter_path)

    # Start the Data Acquisition's thread.
    trigger.execute()
    time.sleep(2*buffer_size)

    # Generate some pulses on the signal outputs by changing the signal output
    # mixer's amplitude. This is for demonstration only and is not necessary to
    # configure the module, we simply generate a signal upon which we can trigger.
    for _ in range(num_pulses):
        daq.setDouble('/%s/sigouts/%d/amplitudes/%d' % (device, out_channel, out_mixer_channel),
                      sigouts_low*(1 + 0.05*np.random.uniform(-1, 1, 1)[0]))
        daq.sync()
        time.sleep(0.2)
        daq.setDouble('/%s/sigouts/%d/amplitudes/%d' % (device, out_channel, out_mixer_channel),
                      sigouts_high*(1 + 0.05*np.random.uniform(-1, 1, 1)[0]))
        daq.sync()
        time.sleep(0.1)
        # Check and display the progress.
        progress = trigger.progress()
        print("Data Acquisition Module progress (acquiring {:d} triggers): {:.2%}.".format(
            trigger_count, progress[0]), end="\r")
        # Check whether the Data Acquisition Module has finished.
        if trigger.finished():
            print("\nTrigger is finished.")
            break
    print("")
    daq.setDouble('/%s/sigouts/%d/amplitudes/%d' % (device, out_channel, out_mixer_channel), sigouts_low)

    # Wait for the Data Acquisition's buffers to finish processing the triggers.
    time.sleep(2*buffer_size)

    # Read the Data Acquisition's data, this command can also be executed before
    # trigger.finished() is True. In that case data recorded up to that point in
    # time is returned and we would still need to issue read() at the end to
    # fetch the rest of the data.
    return_flat_data_dict = True
    data = trigger.read(return_flat_data_dict)
    # Stop the Module's thread and clear the memory.
    trigger.clear()

    # Check that the dictionary returned is non-empty.
    assert data, "read() returned an empty data dictionary, did you subscribe to any paths?"
    # Note: data could be empty if no data arrived, e.g., if the demods were
    # disabled or had rate 0
    assert signal_path in data, "no data recorded: data dictionary has no key `{}`.".format(signal_path)
    samples = data[signal_path]
    filters = data[filter_path]
    print("Data Acquisition's read() returned {} signal segments.".format(len(samples)))
    assert len(samples) == trigger_count, \
        "Unexpected number of signal segments returned: `{}`. Expected: `{}`.".format(len(samples), trigger_count)

    if do_plot and samples:

        import matplotlib.pyplot as plt
        plt.clf()
        # Plot the FFT bins returned by the Data Acquisition.
        for index, sample in enumerate(samples):
            filter = filters[index]
            bin_count = len(sample['value'][0])
            bin_resolution = sample['header']['gridcoldelta']
            frequencies = np.arange(bin_count)
            # Center frequency and bandwidth not yet implemented. So we calculate from the gridcoldelta.
            bandwidth = bin_resolution*len(frequencies)
            frequencies = frequencies*bin_resolution - bandwidth/2.0 + bin_resolution/2.0
            rDb = 20*np.log10(sample['value'][0]*np.sqrt(2)/amplitude)
            rDbCompensated = 20*np.log10((sample['value'][0]/filter['value'][0])*np.sqrt(2)/amplitude)
            plt.subplot(211)
            plt.plot(frequencies, rDb)
            plt.subplot(212)
            plt.plot(frequencies, rDbCompensated)
        plt.subplot(211)
        plt.grid(True)
        title = "Data Acquisition's read() returned {} FFTs ".format(len(samples)) + \
                "each with {} bins".format(len(samples[0]['value'][0]))
        plt.title(title)
        plt.xlabel('Frequency ($Hz$)')
        plt.ylabel('Amplitude R ($dBV$)')
        plt.autoscale(True, 'both', True)
        plt.subplot(212)
        plt.grid(True)
        plt.xlabel('Frequency ($Hz$)')
        plt.ylabel('Amplitude R ($dBV$)(dBV)\n with Demod Filter Compensation')
        plt.autoscale(True, 'both', True)
        plt.draw()
        plt.show()

    return data
