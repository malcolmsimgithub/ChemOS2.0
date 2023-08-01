import json, time

COM_FILE = r'C:\Users\AutoLab\Desktop\PythonLab\pylab\instruments\thorlabs_spectrometer\message.txt'
class ThorlabsCCSHack:
    "Need to run Thorlabs OSA and the script for this hack"

    def __init__(self, com_file = COM_FILE):
        self.com_file = com_file

    def message(self, method, arg):

        # send signal
        message = '\n'.join(["1", method, arg])
        with open(self.com_file, 'w') as f:
            f.write(message)

        # wait
        osa = "1"
        while osa != "0\n":
            time.sleep(0.05)
            with open(self.com_file, 'r') as f:
                osa = f.readline()
                
        # read result
        with open(self.com_file, 'r') as f:
            osa = f.readline()
            result = f.read()
        return result

    def set_integration_time(self, seconds):
        self.message(
            method = 'set_integration_time',
            arg = str(seconds),
        )

    def get_serial(self):
        serial = self.message( method = 'get_serial', arg = '0' )
        return serial

    def single_acquisition(self):
        data = self.message( method = 'single_acquisition', arg = "0" )
        data = data.replace("'", '"')
        data = json.loads(data)
        return data

    def stop(self):
        self.message(method='stop', arg = '0')

    

