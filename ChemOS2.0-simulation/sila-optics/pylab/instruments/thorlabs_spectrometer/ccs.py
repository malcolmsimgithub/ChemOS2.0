from .usb import ModifiedUSBInstrument
from .tlccs_h import *
import numpy as np
import time


## helper functions
TO_STATUS_CODE = {
    'idle': TLCCS_STATUS_SCAN_IDLE,
    'data': TLCCS_STATUS_SCAN_TRANSFER,
    'triggered': TLCCS_STATUS_SCAN_TRIGGERED,
    'start_trans': TLCCS_STATUS_SCAN_START_TRANS,
    'wait_trig': TLCCS_STATUS_WAIT_FOR_EXT_TRIG,
}

def to_status(check):
    return TO_STATUS_CODE.get(check.lower(), TLCCS_STATUS_SCAN_IDLE)

def crc16_block(data):
    crc = 0xFFFF
    data = np.array(data).tobytes()
    for d in data:
        crc = crc16_update(crc, d)
    return crc
def crc16_update(crc, a):
    crc = crc ^ a
    for _ in range(8):
        if crc & 1:
            crc = (crc >> 1) ^ 0xA001
        else:
            crc = (crc >> 1)
    return crc



class ThorlabsCCS(ModifiedUSBInstrument):
    def __init__(self, visa, *args, **kwargs):
        super().__init__(visa, *args, **kwargs)
        self.load_amplitude_correction()
        # self.set_visa_attribute(VI_ATTR_TMO_VALUE, TLCCS_TIMEOUT_DEF)
        self.set_visa_attribute(VI_ATTR_TMO_VALUE, TLCCS_TIMEOUT_MAX)

        # get device parameters
        self.load_wavelength_parameters()
        self.load_dark_current_offset()

        # get device info
        self.model = self.get_visa_attribute(VI_ATTR_MODEL_NAME)
        self.manufacture = self.get_visa_attribute(VI_ATTR_MANF_NAME)
        self.serial = self.get_visa_attribute(VI_ATTR_USB_SERIAL_NUM)
        self.firmware_version, err = self.usb_read(
            TLCCS_RCMD_PRODUCT_INFO, TLCCS_FIRMWARE_VERSION, 0, TLCCS_NUM_VERSION_BYTES * sizeof(ViUInt8), 'u1')
        self.hardware_version, err = self.usb_read(
            TLCCS_RCMD_PRODUCT_INFO, TLCCS_HARDWARE_VERSION, 0, TLCCS_NUM_VERSION_BYTES * sizeof(ViUInt8), 'u1')

        return


    # ---------- intrument functions ------------------------------

    def check_status(self, check = 'idle'):
        status, err = self.get_device_status()

        if isinstance(check, str):
            return bool(status & to_status(check))
        else:
            check = np.array([to_status(c) for c in check], dtype='u2')
            return (status & check).astype(bool)

    def measure(self, repeats = 1, integration_time = None):
        # TODO
        # status, err = self.get_device_status()
        # print(status)

        # only start if it is idle
        # if not self.check_status('idle'):
        #     return np.array([]), StatusCode.error_resource_busy 

        # set integration time
        if integration_time is not None:
            self.set_integration_time(integration_time)

        data = 0.
        self.start_scan_continuous()
        for i in range(repeats):
            # while not self.check_status('data'): # might not need this check
            #     time.sleep(1e-4) # change 1e-5?
            data += self.get_scan_data()[0]
        self.get_integration_time() # use this command to stop
        return data/repeats, StatusCode.success


    def measure_trace(self, repeats = 1, integration_time = None):
        # only start if it is idle
        if not self.check_status('idle'):
            return np.array([]), StatusCode.error_resource_busy 

        # set integration time
        if integration_time is not None:
            self.set_integration_time(integration_time)

        data = 0.
        data_trace = []
        time_list = []
        self.start_scan_continuous()
        for i in range(repeats):
            while not self.check_status('data'): # might not need this check
                time.sleep(1e-4) # change 1e-5?
            new_data = self.get_scan_data()[0]
            data += new_data
            data_trace.append(new_data)
            time_list.append(time.time())
            if i == 0:
                time_start = time_list[i]
            time_list[i] -= time_start
        self.get_integration_time() # use this command to stop
        results = {'time' : time_list, 'average' : data/repeats, 'time_trace' : data_trace}
        return results, StatusCode.success

    def device_info(self):
        device = "{} {}, serial: {}".format(self.manufacture, self.model, self.serial)
        driver = "driver: {:d}.{:d}.{:d}".format(TLCCS_VER_MAJOR, TLCCS_VER_MINOR, TLCCS_VER_SUBMINOR)
        firmware = "firmware: {:d}.{:d}.{:d}/1.{:d}.0".format(*self.firmware_version, self.hardware_version[0])
        return ', '.join([device, driver, firmware])





    # ----------- implmented functions from TLCCS.c ---------------
    def reset(self):
        return self.usb_write(TLCCS_WCMD_RESET, 0, 0, "")

    def get_integration_time(self): 
        "Returns integration time in seconds"
        data, err = self.usb_read(TLCCS_RCMD_INTEGRATION_TIME, 0, 0, TLCCS_NUM_INTEG_CTRL_BYTES * sizeof(ViUInt8), '>u2')
        presc, fill, integ = data & 0x0FFF
        integration_time = (integ-fill+8) * 2. ** presc / 1e6
        return integration_time, err

    def set_integration_time(self, integration_time):
        "Set integration time in seconds"
        # check if integration time is within range
        if integration_time < TLCCS_MIN_INT_TIME or integration_time > TLCCS_MAX_INT_TIME:
            raise ValueError(
                'integration time out of range, must be between {} and {} seconds'.format(TLCCS_MIN_INT_TIME, TLCCS_MAX_INT_TIME))

        integ = int(integration_time * 1e6) # convert to microseconds
        integ_max = int(4095.0 / (1.0 + 0.01 * SH_PERCENT))

        # calculate the prescaler value
        presc = 0
        while integ > integ_max and presc < 20:
            integ = integ >> 1
            presc += 1

        # integ counter has to be so big that CPLD internal counter has finifshed before integ counter becomes zero
        diff = 0
        if integ < (TLCCS_NUM_RAW_PIXELS >> presc):
            diff = (TLCCS_NUM_RAW_PIXELS >> presc) - integ

        # calculate fill
        fill = int( (integ * SH_PERCENT) / 100 + diff )
        integ = integ - 8 + fill

        # when overflow
        # increase the prescaler and decrease integ and fill accordingly
        if integ > TLCCS_NUM_RAW_PIXELS:
            integ = integ >> 1
            integ -= 4
            fill  = fill >> 1
            presc += 1

        # convert to bytes
        data = np.array([presc, fill | 0x1000, integ | 0x2000], dtype='>u2').tobytes()
        return self.usb_write(TLCCS_WCMD_INTEGRATION_TIME, 0, 0, data)

    def start_scan(self):
        err = self.usb_write(TLCCS_WCMD_MODUS, MODUS_INTERN_SINGLE_SHOT, 0)
        return err

    def start_scan_continuous(self):
        err = self.usb_write(TLCCS_WCMD_MODUS, MODUS_INTERN_CONTINUOUS, 0)
        return err

    def get_device_status(self):
        data, err = self.usb_read(TLCCS_RCMD_GET_STATUS, 0, 0, sizeof(ViUInt16), 'u2')
        return data[0], err
    
    def get_scan_data(self): 
        raw = np.frombuffer(self.manager.read_raw(), 'u2').astype('f')
        dark_com = np.mean(raw[DARK_PIXELS_OFFSET:(DARK_PIXELS_OFFSET+NO_DARK_PIXELS)])

        # check if overexposure
        if dark_com > DARK_LEVEL_THRESHOLD_ADC:
            data = np.zeros((TLCCS_NUM_PIXELS, ))
            return data, VI_ERROR_SCAN_DATA_INVALID
            # return data, 'scan data invalid'

        # calculate data
        norm_com = 1.0 / (MAX_ADC_VALUE - dark_com)
        data = (raw[SCAN_PIXELS_OFFSET:SCAN_PIXELS_OFFSET+TLCCS_NUM_PIXELS] - dark_com) * norm_com
        data *= self.fact_acor
        return data, StatusCode.success

    def read_eeprom_no_crc(self, request_value, index, length, dtype = None):
        return self.usb_read(TLCCS_RCMD_READ_EEPROM, request_value, index, length, dtype)
    def read_eeprom(self, request_value, index, length, dtype = None):
        data, err = self.read_eeprom_no_crc(request_value, index, length, dtype)

        if (not err) and request_value >= EE_SW_VERSION:
            csum = crc16_block(data)
            ees, err = self.read_eeprom_no_crc(request_value+length, index, sizeof(ViUInt16), 'u2')
            if csum != ees:
                # print(csum, ees)
                err = VI_ERROR_CYEEPROM_CHKSUM

        return data, err

    def load_amplitude_correction(self): 
        checksums, err = self.read_eeprom(EE_CHECKSUMS, 0, EE_LENGTH_CHECKSUMS, 'u2')

        if err == VI_ERROR_CYEEPROM_CHKSUM:
            chkchk_corrupted_flag = 1
            checksums[TLCCS_ACOR_USER]    = 0
            checksums[TLCCS_ACOR_FACTORY] = 0
            err = StatusCode.success
            print('VI_ERROR_CYEEPROM_CHKSUM encountered')
            # there is something todo here, just incase

        self.fact_acor, err = self.read_eeprom_no_crc(EE_ACOR_FACTORY, 0, EE_LENGTH_ACOR, 'f4')
        self.user_acor, err = self.read_eeprom_no_crc(EE_ACOR_USER, 0, EE_LENGTH_ACOR, 'f4')

        self.fact_sums, err = self.read_eeprom_no_crc(EE_ACOR_FACTORY+EE_LENGTH_ACOR, 0, sizeof(ViUInt16), 'u2')
        self.user_sums, err = self.read_eeprom_no_crc(EE_ACOR_USER+EE_LENGTH_ACOR, 0, sizeof(ViUInt16), 'u2')

        # check if checksums match
        if crc16_block(self.fact_acor) != self.fact_sums:
            return VI_ERROR_CYEEPROM_CHKSUM
        if crc16_block(self.user_acor) != self.user_sums:
            return VI_ERROR_CYEEPROM_CHKSUM
        
        if checksums[TLCCS_ACOR_USER] == self.user_sums:
            self.cal_mode = TLCCS_CAL_MODE_THORLABS
        else:
            self.cal_mode = TLCCS_CAL_MODE_USER
        
        return err


    def load_dark_current_offset(self):
        def _read_offset(offset_type):
            tmp, err = self.read_eeprom(offset_type, 0, EE_LENGTH_OFFSET_MAX, 'u2')
            if err == VI_ERROR_CYEEPROM_CHKSUM:
                tmp = 0xFFFF
                err = StatusCode.success
                # print('VI_ERROR_CYEEPROM_CHKSUM encountered')
            return tmp

        self.even_offset_max = _read_offset(EE_EVEN_OFFSET_MAX)
        self.odd_offset_max = _read_offset(EE_ODD_OFFSET_MAX)
        return StatusCode.success


    def load_wavelength_parameters(self):
        # the poly should be self.fact_poly
        self.fact_poly, err = self.read_eeprom(EE_FACT_CAL_COEF_DATA, 0, EE_LENGTH_FACT_CAL_COEF_DATA, 'f8')
        self.fact_wl, err = self.poly_to_wl_array(self.fact_poly)
        self.user_cal = False

        # user wavelenth calibration
        uc_pixels, uc_wl, uc_cnt, err = self.read_eeuser_points()
        if err == StatusCode.success:
            user_poly = np.polyfit(uc_pixels, uc_wl, 3)[::-1]
            user_wl, err = self.poly_to_wl_array(self.fact_poly)
        if err == StatusCode.success:
            self.user_poly = user_poly
            self.user_wl = user_wl
            self.user_cal = True

        return err

    def get_wavelength_data(self, calibration = 'factory'):
        if calibration.upper() == 'FACTORY':
            return np.copy(self.fact_wl[:TLCCS_NUM_PIXELS]), StatusCode.success
        elif calibration.upper() == 'USER':
            if not self.user_cal:
                return VI_ERROR_TLCCS_NO_USER_DATA
            return np.copy(self.user_wl[:TLCCS_NUM_PIXELS]), StatusCode.success
        else:
            return VI_ERROR_INV_PARAMETER

    def poly_to_wl_array(self, poly):
        pixels = np.arange(0, TLCCS_NUM_PIXELS)
        wl = poly[0] + poly[1] * pixels + poly[2] * pixels ** 2 + poly[3] * pixels ** 3

        # checking if the wavelengths are legit
        if wl[0] < wl[1]:
            if np.sum(wl[1:] - wl[:-1] < 0) != 0: # has wl diff under 0
                return VI_ERROR_TLCCS_INV_USER_DATA
        elif wl[0] > wl[0]:
            if np.sum(wl[1:] - wl[:-1] > 0) != 0: # has wl diff bigger 0
                return VI_ERROR_TLCCS_INV_USER_DATA
        else:
            return VI_ERROR_TLCCS_INV_USER_DATA
        
        return wl, StatusCode.success

    def read_eeuser_points(self):
        cnt, err = self.read_eeprom(EE_USER_CAL_POINTS_CNT, 0, EE_LENGTH_USER_CAL_POINTS_CNT, 'u2')
        if err:
            return [], [], cnt, err
        buf, err = self.read_eeprom(EE_USER_CAL_POINTS_DATA, 0, EE_LENGTH_USER_CAL_POINTS_DATA)
        if err:
            return [], [], cnt, err
        pixels = np.frombuffer(buf[:cnt*sizeof(ViUInt32)], 'u4')
        wl_offset = TLCCS_MAX_NUM_USR_ADJ*sizeof(ViUInt32)
        wl = np.frombuffer(buf[wl_offset:wl_offset+cnt*sizeof(ViReal64)], 'f8')
        return pixels, wl, cnt, err


    # obsolete
    def _measure_slow(self, repeats = 1, integration_time = None):
        # only start if it is idle
        if not self.check_status('idle'):
            return np.array([]), StatusCode.error_resource_busy 

        # set integration time
        if integration_time is not None:
            self.set_integration_time(integration_time)

        # wait until measure finished
        data = 0.
        for i in range(repeats):
            self.start_scan()
            while not self.check_status('idle'):
                time.sleep(0.0001)
            if self.check_status('data'):
                data += self.get_scan_data()[0]
            else:
                return np.array([]), StatusCode.error_resource_busy 
        return data / repeats, StatusCode.success

# static ViStatus tlccs_poly2wlArray(tlccs_wl_cal_t *wl)
# {
#    int      i;
#    ViReal64 d = 0.0;
#    int iDirectionFlag = 0; // 1 means increasing, -1 means decreasing, 0 is an error
#    
#    // check if values are decreasing
#    if(wl->wl[0] < wl->wl[1])
#    {
#       iDirectionFlag = 1;  
#    }
#    else if(wl->wl[0] > wl->wl[1])  
#    {
#       iDirectionFlag = -1; 
#    }
#    else
#       return VI_ERROR_TLCCS_INV_USER_DATA;
#    
#    
#    d = wl->wl[0];
#    for(i = 1; i < TLCCS_NUM_PIXELS; i++)
#    {
#       if(iDirectionFlag == 1) // increasing
#       {
#          if(wl->wl[i] <= d)   return VI_ERROR_TLCCS_INV_USER_DATA; 
#       }
#       else
#       {
#          if(wl->wl[i] >= d)   return VI_ERROR_TLCCS_INV_USER_DATA;
#       }
#
#       d = wl->wl[i];
#    }
#    
#    if(iDirectionFlag == 1)
#    {
#       wl->min     = wl->poly[0];
#       wl->max     = wl->wl[TLCCS_NUM_PIXELS - 1];
#    }
#    else
#    {
#       wl->min     = wl->wl[TLCCS_NUM_PIXELS - 1]; 
#       wl->max     = wl->poly[0];
#    }
#    
#    return VI_SUCCESS;
# }
   
    
# static ViStatus tlccs_getWavelengthParameters (ViSession vi)
# {
#    tlccs_data_t   *data;
#    ViStatus       err = VI_SUCCESS;
#
#    // get private data
#    if((err = viGetAttribute(vi, VI_ATTR_USER_DATA, &data)) != VI_SUCCESS) return err;
#
#    // set the factory calibration valid flag to false
#    data->factory_cal.valid = 0;
#    
#    // read factory adjustment coefficients from EEPROM
#    tlccs_readEEFactoryPoly(vi, data->factory_cal.poly);
#### static ViStatus tlccs_readEEFactoryPoly(ViSession vi, ViReal64 poly[])
####    err = tlccs_readEEPROM(vi, (ViUInt16)EE_FACT_CAL_COEF_DATA, 0, (ViUInt16)EE_LENGTH_FACT_CAL_COEF_DATA, (ViBuf)poly, &cnt);
#
#    tlccs_poly2wlArray(&(data->factory_cal));
#    
#    // read user adjustment nodes from EEPROM and calculate coefficients and wavelength array
#    data->user_cal.valid = 0;
#    err = tlccs_readEEUserPoints(vi, data->user_points.user_cal_node_pixel, data->user_points.user_cal_node_wl, &(data->user_points.user_cal_node_cnt));

   # err = tlccs_readEEPROM(vi, (ViUInt16)EE_USER_CAL_POINTS_CNT, 0, (ViUInt16)EE_LENGTH_USER_CAL_POINTS_CNT, (ViBuf)cnt, &tmp);
   # err = tlccs_readEEPROM(vi, (ViUInt16)EE_USER_CAL_POINTS_DATA, 0, (ViUInt16)EE_LENGTH_USER_CAL_POINTS_DATA, (ViBuf)buf, &tmp);
   # memcpy(pixel, buf, (*cnt) * sizeof(ViUInt32));
   # memcpy(wl, &buf[TLCCS_MAX_NUM_USR_ADJ * sizeof(ViUInt32)], (*cnt) * sizeof(ViReal64));





#### static ViStatus tlccs_nodes2poly(ViInt32 pixel[], ViReal64 wl[], ViInt32 cnt, ViReal64 poly[])
#### {
####    if(LeastSquareInterpolation ((int *)pixel, (double *)wl, (int)cnt, (double *)poly)) return VI_ERROR_TLCCS_INV_USER_DATA;
#    if(err == VI_SUCCESS) err = tlccs_nodes2poly(data->user_points.user_cal_node_pixel, data->user_points.user_cal_node_wl, data->user_points.user_cal_node_cnt, data->user_cal.poly);
#    if(err == VI_SUCCESS) err = tlccs_poly2wlArray(&(data->user_cal));
#    if(err == VI_SUCCESS) data->user_cal.valid = 1;
#                                         
#    return VI_SUCCESS;   // errors ignored by intention
# }

    # def get_amplitude_data(self): return
    # def get_wavelength_data(self): return


    ## thorlabs methods
    # def set_integration_time(self): return
    # def get_integration_time(self): return
    # def start_scan(self): return
    # def start_scan_cont(self): return
    # def start_scan_ext_trg(self): return
    # def start_scan_cont_ext_trg(self): return
    # def get_device_status(self): return
    # def get_scan_data(self): return
    # def get_raw_scan_data(self): return
    # def set_wavelength_data(self): return
    # def get_wavelength_data(self): return
    # def get_user_calibration_points(self): return
    # def set_amplitude_data(self): return
    # def get_amplitude_data(self): return
    # def identification_query(self): return
    # def revision_query(self): return
    # def reset(self): return
    # def self_test(self): return
    # def error_query(self): return
    # def error_message(self): return
    # def set_user_text(self): return
    # def get_user_text(self): return
    # def set_attribute(self): return
    # def get_attribute(self): return
