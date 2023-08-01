from HPLCMS_auto_run_ESI import Auto_HPLCMS
import subprocess
import threading
import time
from run_socket_server import run_socket_server



socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()

server_cmd = f"python -m hplcSilaPackage --ip-address 127.0.0.1 -p 65050 --insecure"

proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)

"""
Default parameters are defined in default_chemspeed.py or defaults.py depending on the mode you are running in.
If you want to change the parameters, you can specify any of those in here.
"""

import datetime
today = datetime.date.today()
today_formatted = today.strftime("%Y%m%d")

params = {
    'Chromeleon_path': {
        'sequence_name': '%s_Autorun_ChemSpeed' %today_formatted,
        'instrument_method_characterization': '/Chemspeed/Instrument_methods/Gradient_Optimized_50%-1.0min_100%-8.0min_100%-14.0min_0.4mL_ESI-pos.instmeth'
        },
    'device_settings': {
        'DADtoLoop_offset': 0.22,
        'DADtoMS_offset': 0.50,
        'splitting_ratio': 0.9,
        'buffer_time': 1.2
        },
    'characterization_setting': {
        'peak_selection': 'biggest',
        },
    'DAD_peak_search_setting': {
        'height_threshold_abs': 2.5,
        'height_threshold_rel': 0.00001,
        'time_range': [0, 15],
        'spectral_range': [250, 500],
        'spectral_range_characterization': [250, 500],
        'spectral_range_for_matching': [250, 500],
        'retentiontime_tolerance': 0.1
    },
    'MS_peak_search_setting_XIC': {
        'low_pass_cutoff': 0.1,
        'mass_tolerance_positive': 0.01,
        'mass_tolerance_negative': 0.01,
        'height_threshold_abs': 1E5,
        'isotope_threshold': None
        }
    }


"""
Initialize the API. if using chemspeed mode, chemspeed should be True and HPLC-MS configuration should be "UofT_LCMS_QE_Chemspeed". If using autosampler mode, chemspeed should be false and configuration should be "UofT_LCMS_QE" or something like that.
"""
LCMS = Auto_HPLCMS(chemspeed=True, **params)
LCMS.auto_run_with_chemspeed(update_retention_list=True)

