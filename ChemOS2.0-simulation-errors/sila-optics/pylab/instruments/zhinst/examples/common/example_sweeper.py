# -*- coding: utf-8 -*-
"""
Zurich Instruments LabOne Python API Example

Demonstrate how to perform a simple frequency sweep using the ziDAQSweeper
class/Sweeper Module.
"""

# Copyright 2016 Zurich Instruments AG

from __future__ import print_function
import time
import numpy as np
#import zhinst.utils
from  pylab.instruments.zhinst import utils


def run_example(device_id, amplitude=0.1, do_plot=False):
    """
    Run the example: Perform a frequency sweep and record demodulator data using
    ziPython's ziDAQSweeper module.

    Requirements:

      Hardware configuration: Connect signal output 1 to signal input 1 with a
      BNC cable.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev2006` or `uhf-dev2006`.

      amplitude (float, optional): The amplitude to set on the signal output.

      do_plot (bool, optional): Specify whether to plot the sweep. Default is no
        plot output.

    Returns:

      sample (list of dict): A list of demodulator sample dictionaries. Each
        entry in the list correspond to the result of a single sweep and is a
        dict containing a demodulator sample.

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
    (daq, device, props) = zhinst.utils.create_api_session(device_id, apilevel_example,
                                                           required_devtype='.*LI|.*IA|.*IS',
                                                           required_err_msg=err_msg)
    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # Now configure the instrument for this experiment. The following channels
    # and indices work on all device configurations. The values below may be
    # changed if the instrument has multiple input/output channels and/or either
    # the Multifrequency or Multidemodulator options installed.
    out_channel = 0
    out_mixer_channel = zhinst.utils.default_output_mixer_channel(props)
    in_channel = 0
    demod_index = 0
    osc_index = 0
    demod_rate = 10e3
    time_constant = 0.01
    exp_setting = [['/%s/sigins/%d/ac'             % (device, in_channel), 0],
                   ['/%s/sigins/%d/range'          % (device, in_channel), 2*amplitude],
                   ['/%s/demods/%d/enable'         % (device, demod_index), 1],
                   ['/%s/demods/%d/rate'           % (device, demod_index), demod_rate],
                   ['/%s/demods/%d/adcselect'      % (device, demod_index), in_channel],
                   ['/%s/demods/%d/order'          % (device, demod_index), 4],
                   ['/%s/demods/%d/timeconstant'   % (device, demod_index), time_constant],
                   ['/%s/demods/%d/oscselect'      % (device, demod_index), osc_index],
                   ['/%s/demods/%d/harmonic'       % (device, demod_index), 1],
                   ['/%s/sigouts/%d/on'            % (device, out_channel), 1],
                   ['/%s/sigouts/%d/enables/%d'    % (device, out_channel, out_mixer_channel), 1],
                   ['/%s/sigouts/%d/range'         % (device, out_channel), 1],
                   ['/%s/sigouts/%d/amplitudes/%d' % (device, out_channel, out_mixer_channel), amplitude]]
    daq.set(exp_setting)

    # Perform a global synchronisation between the device and the data server:
    # Ensure that 1. the settings have taken effect on the device before issuing
    # the poll() command and 2. clear the API's data buffers.
    daq.sync()

    # Create an instance of the Sweeper Module (ziDAQSweeper class).
    sweeper = daq.sweep()

    # Configure the Sweeper Module's parameters.
    # Set the device that will be used for the sweep - this parameter must be set.
    sweeper.set('sweep/device', device)
    # Specify the `gridnode`: The instrument node that we will sweep, the device
    # setting corresponding to this node path will be changed by the sweeper.
    sweeper.set('sweep/gridnode', 'oscs/%d/freq' % osc_index)
    # Set the `start` and `stop` values of the gridnode value interval we will use in the sweep.
    sweeper.set('sweep/start', 4e3)
    if props['devicetype'].startswith('MF'):
        stop = 500e3
    else:
        stop = 50e6
    sweeper.set('sweep/stop', stop)
    # Set the number of points to use for the sweep, the number of gridnode
    # setting values will use in the interval (`start`, `stop`).
    samplecount = 100
    sweeper.set('sweep/samplecount', samplecount)
    # Specify logarithmic spacing for the values in the sweep interval.
    sweeper.set('sweep/xmapping', 1)
    # Automatically control the demodulator bandwidth/time constants used.
    # 0=manual, 1=fixed, 2=auto
    # Note: to use manual and fixed, sweep/bandwidth has to be set to a value > 0.
    sweeper.set('sweep/bandwidthcontrol', 2)
    # Sets the bandwidth overlap mode (default 0). If enabled, the bandwidth of
    # a sweep point may overlap with the frequency of neighboring sweep
    # points. The effective bandwidth is only limited by the maximal bandwidth
    # setting and omega suppression. As a result, the bandwidth is independent
    # of the number of sweep points. For frequency response analysis bandwidth
    # overlap should be enabled to achieve maximal sweep speed (default: 0). 0 =
    # Disable, 1 = Enable.
    sweeper.set('sweep/bandwidthoverlap', 0)

    # Sequential scanning mode (as opposed to binary or bidirectional).
    sweeper.set('sweep/scan', 0)
    # Specify the number of sweeps to perform back-to-back.
    loopcount = 2
    sweeper.set('sweep/loopcount', loopcount)
    # We don't require a fixed sweep/settling/time since there is no DUT
    # involved in this example's setup (only a simple feedback cable), so we set
    # this to zero. We need only wait for the filter response to settle,
    # specified via sweep/settling/inaccuracy.
    sweeper.set('sweep/settling/time', 0)
    # The sweep/settling/inaccuracy' parameter defines the settling time the
    # sweeper should wait before changing a sweep parameter and recording the next
    # sweep data point. The settling time is calculated from the specified
    # proportion of a step response function that should remain. The value
    # provided here, 0.001, is appropriate for fast and reasonably accurate
    # amplitude measurements. For precise noise measurements it should be set to
    # ~100n.
    # Note: The actual time the sweeper waits before recording data is the maximum
    # time specified by sweep/settling/time and defined by
    # sweep/settling/inaccuracy.
    sweeper.set('sweep/settling/inaccuracy', 0.001)
    # Set the minimum time to record and average data to 10 demodulator
    # filter time constants.
    sweeper.set('sweep/averaging/tc', 10)
    # Minimal number of samples that we want to record and average is 100. Note,
    # the number of samples used for averaging will be the maximum number of
    # samples specified by either sweep/averaging/tc or sweep/averaging/sample.
    sweeper.set('sweep/averaging/sample', 10)

    # Now subscribe to the nodes from which data will be recorded. Note, this is
    # not the subscribe from ziDAQServer; it is a Module subscribe. The Sweeper
    # Module needs to subscribe to the nodes it will return data for.x
    path = '/%s/demods/%d/sample' % (device, demod_index)
    sweeper.subscribe(path)

    # Start the Sweeper's thread.
    sweeper.execute()

    start = time.time()
    timeout = 60  # [s]
    print("Will perform", loopcount, "sweeps...")
    while not sweeper.finished():  # Wait until the sweep is complete, with timeout.
        time.sleep(0.2)
        progress = sweeper.progress()
        print("Individual sweep progress: {:.2%}.".format(progress[0]), end="\r")
        # Here we could read intermediate data via:
        # data = sweeper.read(True)...
        # and process it while the sweep is completing.
        # if device in data:
        # ...
        if (time.time() - start) > timeout:
            # If for some reason the sweep is blocking, force the end of the
            # measurement.
            print("\nSweep still not finished, forcing finish...")
            sweeper.finish()
    print("")

    # Read the sweep data. This command can also be executed whilst sweeping
    # (before finished() is True), in this case sweep data up to that time point
    # is returned. It's still necessary still need to issue read() at the end to
    # fetch the rest.
    return_flat_dict = True
    data = sweeper.read(return_flat_dict)
    sweeper.unsubscribe(path)

    # Stop the sweeper thread and clear the memory.
    sweeper.clear()

    # Check the dictionary returned is non-empty.
    assert data, "read() returned an empty data dictionary, did you subscribe to any paths?"
    # Note: data could be empty if no data arrived, e.g., if the demods were
    # disabled or had rate 0.
    assert path in data, "No sweep data in data dictionary: it has no key '%s'" % path
    samples = data[path]
    print("Returned sweeper data contains", len(samples), "sweeps.")
    assert len(samples) == loopcount, \
        "The sweeper returned an unexpected number of sweeps: `%d`. Expected: `%d`." % (len(samples), loopcount)

    if do_plot:
        import matplotlib.pyplot as plt
        _, (ax1, ax2) = plt.subplots(2, 1)

        for sample in samples:
            frequency = sample[0]['frequency']
            R = np.abs(sample[0]['x'] + 1j*sample[0]['y'])
            phi = np.angle(sample[0]['x'] + 1j*sample[0]['y'])
            ax1.plot(frequency, R)
            ax2.plot(frequency, phi)
        ax1.set_title('Results of %d sweeps.' % len(samples))
        ax1.grid()
        ax1.set_ylabel(r'Demodulator R ($V_\mathrm{RMS}$)')
        ax1.set_xscale('log')
        ax1.set_ylim(0.0, 0.1)

        ax2.grid()
        ax2.set_xlabel('Frequency ($Hz$)')
        ax2.set_ylabel(r'Demodulator Phi (radians)')
        ax2.set_xscale('log')
        ax2.autoscale()

        plt.draw()
        plt.show()

    return samples
