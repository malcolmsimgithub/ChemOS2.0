import clr, time
import ThorlabsSpectrum as TS
COM_FILE = r'C:\Users\AutoLab\Desktop\PythonLab\pylab\instruments\thorlabs_spectrometer\message.txt'

class CCS:
    def __init__(self):
        self.ccs = _deviceMgr.GetDevice[TS.CCS](0)
        
    def get_serial(self, *args):
        return self.ccs.GetSerial()
        
    def set_integration_time(self, seconds, *args):
        self.ccs.IntegrationTime_ms.Set(seconds*1000)
        return ''
        
    def single_acquisition(self, *args):
        # acquisition
        self.ccs.StartSingleAcquisition()
        while self.ccs.IsRunningAcquisition():
            time.sleep(0.02)

        # getspectrum
        spectrum = TS.Data.TraceData_Spectrum()
        self.ccs.GetSpectrum(spectrum)
        x = list(spectrum.GetSpectrum().X)
        y = list(spectrum.GetSpectrum().Value)
        return {'x': x, 'y': y}
        


class Runner:
    def __init__(self, com_file, obj):
        self.com_file = com_file
        self.obj = obj

        with open(self.com_file, 'w') as f:
            f.write('0\n1\n1')
    
    def wait_message(self):

        # wait
        osa = "0"
        while osa != "1\n":
            time.sleep(0.1)
            with open(self.com_file, 'r') as f:
                osa = f.readline()

        # read message
        with open(self.com_file, 'r') as f:
            osa = f.readline()
            method = f.readline().rstrip()
            arg = f.readline().rstrip()
        
        # stop
        if method == 'stop':
            with open(self.com_file, 'w') as f:
                f.write('0\n{}'.format(""))
            return False
        
        # Execute
        print('Execute {}({})'.format(method, arg))
        result = getattr(self.obj, method)(float(arg))
        with open(self.com_file, 'w') as f:
            f.write('0\n{}'.format(str(result)))
        return True
    
    def listen(self):
        loop = True
        while loop:
            loop = self.wait_message()
            time.sleep(0.05)


        


# Run
ccs = CCS()
run = Runner(com_file = COM_FILE, obj = ccs)
run.listen()

# self.ccs.ResetDevice
# self.ccs.GetInstrumentHealth
