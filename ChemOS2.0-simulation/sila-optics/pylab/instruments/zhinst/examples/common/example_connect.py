# -*- coding: utf-8 -*-
"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments device via the Data Server
program.
"""

# Copyright 2016 Zurich Instruments AG

from __future__ import print_function
import numpy as np
#import zhinst.utils
from  pylab.instruments.zhinst import utils

def run_example(device_id):
    """Run the example: Create an API session by connecting to a Zurich Instruments
    device via the Data Server, ensure the demodulators are enabled and obtain a
    single demodulator sample via getSample(). Calculate the sample's RMS
    amplitude and add it as a field to the "sample" dictionary.

    Note:

      This is intended to be a simple example demonstrating how to connect to a
      Zurich Instruments device from ziPython. In most cases, data acquisition
      should use either ziDAQServer's poll() method or an instance of the
      ziDAQRecorder class, not the getSample() method.

    Requirements:

      Hardware configuration: Connect signal output 1 to signal input 1 with a
      BNC cable.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev2006` or `uhf-dev2006`.

    Returns:

      sample (dict): The acquired demodulator sample dictionary.

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
    err_msg = "This example only supports instruments with demodulators."
    (daq, device, _) = zhinst.utils.create_api_session(device_id, apilevel_example,
                                                       required_devtype='.*LI|.*IA|.*IS',
                                                       required_err_msg=err_msg)
    zhinst.utils.api_server_version_check(daq)

    # Enable the demodulator.
    daq.setInt('/%s/demods/0/enable' % device, 1)
    # Set the demodulator output rate.
    daq.setDouble('/%s/demods/0/rate' % device, 1.0e3)

    # Perform a global synchronisation between the device and the data server:
    # Ensure that the settings have taken effect on the device before issuing
    # the getSample() command.
    daq.sync()

    # Obtain one demodulator sample. If the demodulator is not enabled (as
    # above) then the command will time out: we'll get a RuntimeError showing
    # that a `ZIAPITimeoutException` occured.
    sample = daq.getSample('/%s/demods/0/sample' % device)
    # Calculate the demodulator's magnitude and phase and add them to the sample
    # dict.
    sample['R'] = np.abs(sample['x'] + 1j*sample['y'])
    sample['phi'] = np.angle(sample['x'] + 1j*sample['y'])
    print("Measured RMS amplitude is {:.3e} V.".format(sample['R'][0]))
    return sample
