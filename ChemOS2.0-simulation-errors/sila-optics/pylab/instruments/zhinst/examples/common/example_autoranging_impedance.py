# -*- coding: utf-8 -*-
"""
Zurich Instruments LabOne Python API Example

Demonstrate how to perform a manually triggered autoranging for impedance while working in
manual range mode.
"""

# Copyright 2017 Zurich Instruments AG

from __future__ import print_function
import time
#import zhinst.utils
from  pylab.instruments.zhinst import utils

def run_example(device_id):
    """
    Run the example: Set the device to impedance manual range mode and execute
    a single auto ranging event.

    Requirements:

      Hardware configuration: Connect Impedance Testfixture and attach a 1kOhm
      resistor

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev3300` or `mf-dev3300`.

    Returns:

      true

    Raises:

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programming Manual" for further help, available:
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
    err_msg = "This example only supports instruments with IA option."
    (daq, device, _) = zhinst.utils.create_api_session(device_id, apilevel_example,
                                                       required_options=['IA'],
                                                       required_err_msg=err_msg)
    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # Now configure the instrument for this experiment. The following channels
    # and indices work on all devices with IA option. The values below may be
    # changed if the instrument has multiple IA modules.
    imp_index = 0
    curr_index = daq.getInt('/%s/imps/%d/current/inputselect' % (device, imp_index))
    volt_index = daq.getInt('/%s/imps/%d/voltage/inputselect' % (device, imp_index))
    man_curr_range = 10e-3
    man_volt_range = 10e-3
    exp_settings = [['/%s/imps/%d/enable' % (device, imp_index), 1],
                    ['/%s/imps/%d/mode' % (device, imp_index), 0],
                    ['/%s/imps/%d/auto/output' % (device, imp_index), 1],
                    ['/%s/imps/%d/auto/bw' % (device, imp_index), 1],
                    ['/%s/imps/%d/freq' % (device, imp_index), 500],
                    ['/%s/imps/%d/auto/inputrange' % (device, imp_index), 0],
                    ['/%s/currins/%d/range' % (device, curr_index), man_curr_range],
                    ['/%s/sigins/%d/range' % (device, volt_index), man_volt_range]]
    daq.set(exp_settings)
    # Perform a global synchronisation between the device and the data server:
    # Ensure that the settings have taken effect on the device before setting
    # the next configuration.
    daq.sync()

    # After setting the device in manual ranging mode we want to trigger manually
    # a one time auto ranging to find a suitable range. Therefore, we trigger the
    # auto ranging for the current input as well as for the voltage input.
    trigger_auto_ranging = [['/%s/currins/%d/autorange' % (device, curr_index), 1],
                            ['/%s/sigins/%d/autorange' % (device, volt_index), 1]]
    print('Start auto ranging. This takes a few seconds.')
    daq.set(trigger_auto_ranging)

    t_start = time.time()
    timeout = 20
    finished = 0
    # The auto ranging takes some time. We do not want to continue before the
    # best range is found. Therefore, we implement a loop to check if the auto
    # ranging is finished. These nodes maintain value 1 until autoranging has
    # finished.
    while not finished:
        time.sleep(0.5)
        currins_autorange = daq.getInt('/%s/currins/%d/autorange' % (device, curr_index))
        sigins_autorange = daq.getInt('/%s/sigins/%d/autorange'  % (device, volt_index))
        # We are finished when both nodes have been set back to 0 by the device.
        finished = (currins_autorange == 0) and (sigins_autorange == 0)
        if time.time() - t_start > timeout:
            raise Exception("Autoranging failed after {} seconds.".format(timeout))
    print('Auto ranging finished after {:0.1f} s.'.format(time.time() - t_start))

    auto_curr_range = daq.getDouble('/%s/currins/%d/range' % (device, curr_index))
    auto_volt_range = daq.getDouble('/%s/sigins/%d/range' % (device, volt_index))
    print('Current range changed from {:0.1e} A to {:0.1e} A.'.format(man_curr_range, auto_curr_range))
    print('Voltage range changed from {:0.1e} V to {:0.1e} V.'.format(man_volt_range, auto_volt_range))

    return True
