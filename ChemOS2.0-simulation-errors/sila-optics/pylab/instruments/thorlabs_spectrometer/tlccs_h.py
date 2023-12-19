from pyvisa.constants import *
from pyvisa.ctwrapper.types import *
from ctypes import sizeof

# Mine
_VI_ERROR                   = 0x80000000

# FROM TLCCS.h
#ifndef __TLCCS_H__
__TLCCS_H__                      = 00                              
#if defined(__cplusplus) || defined(__cplusplus__)
#endif
#ifndef ViAttrState
#if defined(_VI_INT64_UINT64_DEFINED) && defined(_VISA_ENV_IS_64_BIT)
#else
#endif
#endif
TLCCS_FIND_PATTERN               = "USB?*?{VI_ATTR_MANF_ID==0x1313 && ((VI_ATTR_MODEL_CODE==0x8081) || (VI_ATTR_MODEL_CODE==0x8083) || (VI_ATTR_MODEL_CODE==0x8085) || (VI_ATTR_MODEL_CODE==0x8087) || (VI_ATTR_MODEL_CODE==0x8089))}"
TLCCS_VID                        = 0x1313                            # Thorlabs
CCS100_PID                       = 0x8081                            # CCS100 Compact Spectrometer
CCS125_PID                       = 0x8083                            # CCS125 Special Spectrometer
CCS150_PID                       = 0x8085                            # CCS150 UV Spectrometer
CCS175_PID                       = 0x8087                            # CCS175 NIR Spectrometer
CCS200_PID                       = 0x8089                            # CCS200 UV-NIR Spectrometer
TLCCS_EXTRACT_MAJOR              = lambda revision:                ((revision & 0xFFF00000) >> 20)
TLCCS_EXTRACT_MINOR              = lambda revision:                ((revision & 0x000FFF00) >> 8)
TLCCS_EXTRACT_SUBMINOR           = lambda revision:             (revision & 0x000000FF)
TLCCS_TIMEOUT_MIN                = 1000                            
TLCCS_TIMEOUT_DEF                = 2000                            
TLCCS_TIMEOUT_MAX                = 60000                           
TLCCS_BUFFER_SIZE                = 256                               # general buffer size
TLCCS_ERR_DESCR_BUFFER_SIZE      = 512                               # buffer size for error messages
TLCCS_TEXT_BUFFER_SIZE           = TLCCS_BUFFER_SIZE                 # buffer size for texts from the SPX
TLCCS_NUM_PIXELS                 = 3648                              # number of effective pixels of CCD
TLCCS_NUM_RAW_PIXELS             = 3694                              # number of raw pixels
TLCCS_MAX_USER_NAME_SIZE         = 32                                # including the trailing '\0'
TLCCS_MIN_NUM_USR_ADJ            = 4                                 # minimum number of user adjustment data points
TLCCS_MAX_NUM_USR_ADJ            = 10                                # maximum number of user adjustment data points
TLCCS_ATTR_USER_DATA             = (0x3FFF0007)                      # we do not recursively resolve this value because the LabView driver importer cannot either
TLCCS_ATTR_CAL_MODE              = (0x3FFA0000)                      # this attribute is either 0 = user or 1 = THORLABS amplitude correction
TLCCS_MAX_INT_TIME               = 60.0                              # 60s is the maximum integration time
TLCCS_MIN_INT_TIME               = 0.00001                           # 10us is the minimum integration time
TLCCS_DEF_INT_TIME               = 0.01                              # 10ms is the default integration time
TLCCS_CAL_MODE_USER              = 0                               
TLCCS_CAL_MODE_THORLABS          = 1                               
TLCCS_AMP_CORR_FACT_MIN          = 0.001                             # the minimum correction factor
TLCCS_AMP_CORR_FACT_MAX          = 1000.0                            # the maximum correction factor, standard is 1.0
VI_ERROR_NSUP_COMMAND            = (_VI_ERROR + 0x3FFC0801)          # 0xBFFC0801
VI_ERROR_TLCCS_UNKNOWN           = (_VI_ERROR + 0x3FFC0802)          # 0xBFFC0802
VI_ERROR_SCAN_DATA_INVALID       = (_VI_ERROR + 0x3FFC0803)          # 0xBFFC0803
VI_ERROR_XSVF_SIZE               = (_VI_ERROR + 0x3FFC0A00)          # 0xBFFC0A00
VI_ERROR_XSVF_MEMORY             = (_VI_ERROR + 0x3FFC0A01)          # 0xBFFC0A01
VI_ERROR_XSVF_FILE               = (_VI_ERROR + 0x3FFC0A02)          # 0xBFFC0A02
VI_ERROR_FIRMWARE_SIZE           = (_VI_ERROR + 0x3FFC0A10)          # 0xBFFC0A10
VI_ERROR_FIRMWARE_MEMORY         = (_VI_ERROR + 0x3FFC0A11)          # 0xBFFC0A11
VI_ERROR_FIRMWARE_FILE           = (_VI_ERROR + 0x3FFC0A12)          # 0xBFFC0A12
VI_ERROR_FIRMWARE_CHKSUM         = (_VI_ERROR + 0x3FFC0A13)          # 0xBFFC0A13
VI_ERROR_FIRMWARE_BUFOFL         = (_VI_ERROR + 0x3FFC0A14)          # 0xBFFC0A14
VI_ERROR_CYEEPROM_SIZE           = (_VI_ERROR + 0x3FFC0A20)          # 0xBFFC0A20
VI_ERROR_CYEEPROM_MEMORY         = (_VI_ERROR + 0x3FFC0A21)          # 0xBFFC0A21
VI_ERROR_CYEEPROM_FILE           = (_VI_ERROR + 0x3FFC0A22)          # 0xBFFC0A22
VI_ERROR_CYEEPROM_CHKSUM         = (_VI_ERROR + 0x3FFC0A23)          # 0xBFFC0A23
VI_ERROR_CYEEPROM_BUFOVL         = (_VI_ERROR + 0x3FFC0A24)          # 0xBFFC0A24
VI_ERROR_USBCOMM_OFFSET          = (_VI_ERROR + 0x3FFC0B00)          # 0xBFFC0B00
VI_ERROR_TLCCS_ENDP0_SIZE        = (VI_ERROR_USBCOMM_OFFSET + 0x01)  # 0xBFFC0B01
VI_ERROR_TLCCS_EEPROM_ADR_TO_BIG = (VI_ERROR_USBCOMM_OFFSET + 0x02)  # 0xBFFC0B02
VI_ERROR_TLCCS_XSVF_UNKNOWN      = (VI_ERROR_USBCOMM_OFFSET + 0x11)  # 0xBFFC0B11
VI_ERROR_TLCCS_XSVF_TDOMISMATCH  = (VI_ERROR_USBCOMM_OFFSET + 0x12)  # 0xBFFC0B12
VI_ERROR_TLCCS_XSVF_MAXRETRIES   = (VI_ERROR_USBCOMM_OFFSET + 0x13)  # 0xBFFC0B13
VI_ERROR_TLCCS_XSVF_ILLEGALCMD   = (VI_ERROR_USBCOMM_OFFSET + 0x14)  # 0xBFFC0B14
VI_ERROR_TLCCS_XSVF_ILLEGALSTATE = (VI_ERROR_USBCOMM_OFFSET + 0x15)  # 0xBFFC0B15
VI_ERROR_TLCCS_XSVF_DATAOVERFLOW = (VI_ERROR_USBCOMM_OFFSET + 0x16)  # 0xBFFC0B16
VI_ERROR_TLCCS_I2C_NACK          = (VI_ERROR_USBCOMM_OFFSET + 0x20)  # 0xBFFC0B20
VI_ERROR_TLCCS_I2C_ERR           = (VI_ERROR_USBCOMM_OFFSET + 0x21)  # 0xBFFC0B21
VI_ERROR_TLCCS_READ_INCOMPLETE   = (VI_ERROR_USBCOMM_OFFSET + 0x40)
VI_ERROR_TLCCS_NO_USER_DATA      = (VI_ERROR_USBCOMM_OFFSET + 0x41)  # there is no wavelength adjustment data available at the instruments nonvolatile memory
VI_ERROR_TLCCS_INV_USER_DATA     = (VI_ERROR_USBCOMM_OFFSET + 0x42)  # the given wavelength adjustment data results in negative wavelength values.
VI_ERROR_TLCCS_INV_ADJ_DATA      = (VI_ERROR_USBCOMM_OFFSET + 0x43)  # read out amplitude/wavelength adjustment data is out of range/corrupt
TLCCS_STATUS_SCAN_IDLE           = 0x0002                            # CCS waits for new scan to execute
TLCCS_STATUS_SCAN_TRIGGERED      = 0x0004                            # scan in progress
TLCCS_STATUS_SCAN_START_TRANS    = 0x0008                            # scan starting
TLCCS_STATUS_SCAN_TRANSFER       = 0x0010                            # scan is done, waiting for data transfer to PC
TLCCS_STATUS_WAIT_FOR_EXT_TRIG   = 0x0080                            # same as IDLE except that external trigger is armed
TLCCS_CAL_DATA_SET_FACTORY       = 0                               
TLCCS_CAL_DATA_SET_USER          = 1                               
ACOR_APPLY_TO_MEAS               = 1                                 # macro for parameter mode
ACOR_APPLY_TO_MEAS_NVMEM         = 2                                 # macro for parameter mode
ACOR_FROM_CURRENT                = 1                                 # macro for parameter mode
ACOR_FROM_NVMEM                  = 2                                 # macro for parameter mode
#if defined(__cplusplus) || defined(__cplusplus__)
#endif
#endif  /* __TLCCS_H__ */

# FROM TLCCS.c
#ifdef _CVI_
#else
#endif
NIVISA_USB                     = 00                                                             # this is necessary to access NI-VISA-USB
TLCCS_VER_MAJOR                = 2                                                            
TLCCS_VER_MINOR                = 0                                                            
TLCCS_VER_SUBMINOR             = 0                                                            
TLCCS_MAKE_REVISION            = lambda major, minor, subminor:  (((major & 0x00000FFF) << 20) | ((minor & 0x00000FFF) << 8) | (subminor & 0x000000FF))
TLCCS_VERSION                  = TLCCS_MAKE_REVISION(TLCCS_VER_MAJOR, TLCCS_VER_MINOR, TLCCS_VER_SUBMINOR)
# #ifdef _CVI_DEBUG_
# TLCCS_LOCK_STATE               = VI_NULL                                                      
# #else
# TLCCS_LOCK_STATE               = VI_EXCLUSIVE_LOCK                                            
# #endif
TLCCS_SERIAL_NO_LENGTH         = 24                                                           
TLCCS_NUM_POLY_POINTS          = 4                                                            
TLCCS_NUM_INTEG_CTRL_BYTES     = 6                                                            
TLCCS_NUM_VERSION_BYTES        = 3                                                            
TLCCS_NUM_FLAG_WORDS           = 1                                                            
TLCCS_NUM_CHECKSUMS            = 2                                                            
TLCCS_FIRMWARE_VERSION         = 0                                                            
TLCCS_HARDWARE_VERSION         = 1                                                            
TLCCS_ACOR_FACTORY             = 0                                                            
TLCCS_ACOR_USER                = 1                                                            
# TLCCS_DEFAULT_ERR_QUERY_MODE   = VI_ON                                                        
# INVAL_RANGE                    = (val, min, max)     ( ((val) < (min)) || ((val) > (max)) )   
DEBUG_BUF_SIZE                 = 512                                                          
ENDPOINT_0_TRANSFERSIZE        = 64                                                             # this is the max. size of bytes that can be transferred at once for Endpoint 0
MATRIX_ROWS                    = 4                                                            
MATRIX_COLS                    = 4                                                            
EE_LENGTH_SERIAL_NO            = ( sizeof(ViChar)                       * TLCCS_SERIAL_NO_LENGTH    )
EE_LENGTH_SW_VERSION           = (4                                                                 )
EE_LENGTH_USER_LABEL           = ( sizeof(ViChar)                       * TLCCS_MAX_USER_NAME_SIZE  )
EE_LENGTH_FACT_CAL_COEF_FLAG   = (2                                                                 )
EE_LENGTH_FACT_CAL_COEF_DATA   = ( sizeof(ViReal64)                     * TLCCS_NUM_POLY_POINTS     )
EE_LENGTH_USER_CAL_COEF_FLAG   = (2                                                                 )
EE_LENGTH_USER_CAL_COEF_DATA   = ( sizeof(ViReal64)                     * TLCCS_NUM_POLY_POINTS     )
EE_LENGTH_USER_CAL_POINTS_CNT  = (2                                                                 )
EE_LENGTH_USER_CAL_POINTS_DATA = ((sizeof(ViInt32) + sizeof(ViReal64))  * TLCCS_MAX_NUM_USR_ADJ     )
EE_LENGTH_OFFSET_MAX           = (2                                                                 )
EE_LENGTH_ACOR                 = ( sizeof(ViReal32)                     * TLCCS_NUM_PIXELS          )
EE_LENGTH_FLAGS                = ( sizeof(ViUInt32)                     * TLCCS_NUM_FLAG_WORDS      )
EE_LENGTH_CHECKSUMS            = ( sizeof(ViUInt16)                     * TLCCS_NUM_CHECKSUMS       )
EE_SIZE_CHECKSUM               = (2)                                                          
EE_SIZE_BOOT_CODE              = (1)                                                          
EE_SIZE_SERIAL_NO              = (EE_LENGTH_SERIAL_NO)                                        
EE_SIZE_SW_VERSION             = (EE_LENGTH_SW_VERSION            + EE_SIZE_CHECKSUM)         
EE_SIZE_USER_LABEL             = (EE_LENGTH_USER_LABEL            + EE_SIZE_CHECKSUM)         
EE_SIZE_FACT_CAL_COEF_FLAG     = (EE_LENGTH_FACT_CAL_COEF_FLAG    + EE_SIZE_CHECKSUM)         
EE_SIZE_FACT_CAL_COEF_DATA     = (EE_LENGTH_FACT_CAL_COEF_DATA    + EE_SIZE_CHECKSUM)         
EE_SIZE_USER_CAL_COEF_FLAG     = (EE_LENGTH_USER_CAL_COEF_FLAG    + EE_SIZE_CHECKSUM)         
EE_SIZE_USER_CAL_COEF_DATA     = (EE_LENGTH_USER_CAL_COEF_DATA    + EE_SIZE_CHECKSUM)         
EE_SIZE_USER_CAL_POINTS_CNT    = (EE_LENGTH_USER_CAL_POINTS_CNT   + EE_SIZE_CHECKSUM)         
EE_SIZE_USER_CAL_POINTS_DATA   = (EE_LENGTH_USER_CAL_POINTS_DATA  + EE_SIZE_CHECKSUM)         
EE_SIZE_OFFSET_MAX             = (EE_LENGTH_OFFSET_MAX            + EE_SIZE_CHECKSUM)         
EE_SIZE_ACOR                   = (EE_LENGTH_ACOR                  + EE_SIZE_CHECKSUM)         
EE_SIZE_FLAGS                  = (EE_LENGTH_FLAGS                 + EE_SIZE_CHECKSUM)         
EE_SIZE_CHECKSUMS              = (EE_LENGTH_CHECKSUMS             + EE_SIZE_CHECKSUM)         
EE_BOOT_CODE                   = 0                                                            
EE_VENDOR_ID                   = 1                                                            
EE_PRODUCT_ID                  = 3                                                            
EE_DEVICE_ID                   = 5                                                            
EE_SERIAL_NO                   = 8                                                            
EE_SW_VERSION                  = (EE_SERIAL_NO              + EE_SIZE_SERIAL_NO              )  # software version
EE_USER_LABEL                  = (EE_SW_VERSION             + EE_SIZE_SW_VERSION             )  # user label
EE_FACT_CAL_COEF_FLAG          = (EE_USER_LABEL             + EE_SIZE_USER_LABEL             )  # factory calibration flags
EE_FACT_CAL_COEF_DATA          = (EE_FACT_CAL_COEF_FLAG     + EE_SIZE_FACT_CAL_COEF_FLAG     )  # factory calibration coefficients
EE_USER_CAL_COEF_FLAG          = (EE_FACT_CAL_COEF_DATA     + EE_SIZE_FACT_CAL_COEF_DATA     )  # user calibration flags
EE_USER_CAL_COEF_DATA          = (EE_USER_CAL_COEF_FLAG     + EE_SIZE_USER_CAL_COEF_FLAG     )  # user calibration coefficients
EE_USER_CAL_POINTS_CNT         = (EE_USER_CAL_COEF_DATA     + EE_SIZE_USER_CAL_COEF_DATA     )  # user calibration points count
EE_USER_CAL_POINTS_DATA        = (EE_USER_CAL_POINTS_CNT    + EE_SIZE_USER_CAL_POINTS_CNT    )  # user calibration points
EE_EVEN_OFFSET_MAX             = (EE_USER_CAL_POINTS_DATA   + EE_SIZE_USER_CAL_POINTS_DATA   )  # even offset max
EE_ODD_OFFSET_MAX              = (EE_EVEN_OFFSET_MAX        + EE_SIZE_OFFSET_MAX             )  # odd offset max
EE_ACOR_FACTORY                = (EE_ODD_OFFSET_MAX         + EE_SIZE_OFFSET_MAX             )  # amplitude correction, factory setting
EE_ACOR_USER                   = (EE_ACOR_FACTORY           + EE_SIZE_ACOR                   )  # amplitude correction, factory setting
EE_FLAGS                       = (EE_ACOR_USER              + EE_SIZE_ACOR                   )  # flags for e.g. user cal/factory cal
EE_CHECKSUMS                   = (EE_FLAGS                  + EE_SIZE_FLAGS                  )  # checksums for amplitude correction arrays
EE_FREE                        = (EE_CHECKSUMS              + EE_SIZE_CHECKSUMS              )  # free memory
TLCCS_ERR_SRC_EEPROM           = 1                                                            
TLCCS_ERR_SRC_FPGA             = 2                                                            
TLCCS_ERR_SRC_FIRMWARE         = 3                                                            
TLCCS_WCMD_WRITE_EEPROM        = 0x21                                                         
TLCCS_WCMD_INTEGRATION_TIME    = 0x23                                                         
TLCCS_WCMD_MODUS               = 0x24                                                         
TLCCS_WCMD_RESET               = 0x26                                                         
TLCCS_RCMD_READ_EEPROM         = 0x21                                                         
TLCCS_RCMD_INTEGRATION_TIME    = 0x23                                                         
TLCCS_RCMD_PRODUCT_INFO        = 0x25                                                         
TLCCS_RCMD_GET_STATUS          = 0x30                                                         
TLCCS_RCMD_GET_ERROR           = 0xFF                                                         
MODUS_INTERN_SINGLE_SHOT       = 0                                                            
MODUS_INTERN_CONTINUOUS        = 1                                                            
MODUS_EXTERN_SINGLE_SHOT       = 2                                                            
MODUS_EXTERN_CONTINUOUS        = 3                                                            
TLCCS_CALIB_VALID_FLAG         = 0x5A                                                           # this is the value for check bytes
TLCCS_USERCAL_VALID_FLAG       = 0x5A                                                           # user wavelenght adjustment data is valid
DEFAULT_USER_TEXT              = "My CCS Spectrometer"                                        
#ifndef DO_NOT_COMPILE_DRIVER_HERE
#else
#endif
#ifndef DO_NOT_COMPILE_DRIVER_HERE
#else
#endif
#ifndef DO_NOT_COMPILE_DRIVER_HERE
SH_PERCENT                     = 16.5                                                         
MIN_SHUTTER_PULSES             = ((int)(3695 * SH_PERCENT))                           
FACTORY_ADJ_OFFSET             = 62749006                                                     
FACTORY_SET_START              = 19901201                                                     
THORLABS_SET_START             = 91901201                                                     
TARGET_FACTROY_DATA            = 0                                                            
TARGET_USER_DATA               = 1                                                            
TARGET_THORLABS_DATA           = 2                                                            
#ifdef NOT_NEEDED_UNTIL_NOW
#endif   // NOT_NEEDED_UNTIL_NOW
#ifdef RESERVED_FOR_LATER_USE
#endif   // RESERVED_FOR_LATER_USE
NO_DARK_PIXELS                 = 12                                                             # we got 12 dark pixels
DARK_PIXELS_OFFSET             = 16                                                             # dark pixels start at positon 16 within raw data
SCAN_PIXELS_OFFSET             = 32                                                             # real measurement start at position 32 within raw data
MAX_ADC_VALUE                  = 0xFFFF                                                         # this is full scale of a 16bit Analog Digital Converter (ADC)
DARK_LEVEL_THRESHOLD           = (0.99)                                                         # when dark level is above 99% f.s. of ADC mark scan as invalid (overexposed)
DARK_LEVEL_THRESHOLD_ADC       = (DARK_LEVEL_THRESHOLD * (MAX_ADC_VALUE))             
#ifdef NOT_NEEDED_UNTIL_NOW
#endif   // NOT_NEEDED_UNTIL_NOW
#ifdef NOT_NEEDED_UNTIL_NOW
#endif   // NOT_NEEDED_UNTIL_NOW
#ifdef NOT_NEEDED_UNTIL_NOW
#endif   // NOT_NEEDED_UNTIL_NOW
#endif   // DO_NOT_COMPILE_DRIVER_HERE

