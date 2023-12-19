# import base
from .base import *

from .thorlabs_spectrometer.ccs import ThorlabsCCS
from .thorlabs_led_driver import ThorlabsDC4100
from .thorlabs_powermeter import ThorlabsPM100D
from .thorlabs_kinesis import ThorlabsKSC101, ThorlabsK10CR1
from .thorlabs_filter import ThorlabsFW212C

from .pumps.hamilton_pump import PumpPSD8
from .pumps.tecan_pump import XCPump

from .heidolph_hotplate import HotplateHei

from .nanalysis_nmr import NMR60Pro

from .picoquant.tcspc import TH260

from .numato_usbgpio import Numato_Usbgpio
from .numato_usbgpio import Numato_Usbrelay

from .picoscope.ps5242D import PS5242D

from .newport_monochromator import CS210

# from .zurich_lockin import MFLI500

from .Instek_FG import Instek_FG

from .ftdi import TTL232R

from .valco_valve import Valco_valve

from .DH_mini import DH_mini

from .sartorius_balance import EntrisBalance

from .coherent_obis import Coherent_OBIS

from .harverd_pump import Harverd_11Elite

from .Torrey_pains_shaker import TorreyPinesSC20

#from .chemspeed_controller.controller import ChemspeedController
#from .pypowerusb.powerUSB import powerUSB
#from .thorlabs_APT.thorlabs_K10CR1 import ThorlabsK10CR1
