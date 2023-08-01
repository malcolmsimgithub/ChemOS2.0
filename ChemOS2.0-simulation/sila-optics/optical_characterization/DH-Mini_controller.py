import pylab.instruments as ins
from pylab.instruments import pv_const
import numpy as np
import matplotlib.pyplot as plt
import time

#initialize the device
DM= ins.DH_mini('ASRL10::INSTR', verbose=False)

DM.shutter_open()
