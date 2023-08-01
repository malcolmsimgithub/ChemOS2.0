try:
    from .numato_usbgpio import Numato_Usbgpio
except:
    from numato_usbgpio import Numato_Usbgpio

class DH_mini():
    
    def __init__(self, visa, verbose = False):
        self.gpio = Numato_Usbgpio(visa)
        self.verbose = verbose

    def deuterium_on(self):
        self.gpio.gpio_on(5)
        if self.verbose:
            print('#### DH_mini : deuterium lamp ON #############')
    
    def deuterium_off(self):
        self.gpio.gpio_off(5)
        if self.verbose:
            print('#### DH_mini : deuterium lamp OFF #############')

    def halogen_on(self):
        self.gpio.gpio_on(7)
        if self.verbose:
            print('#### DH_mini : halogen lamp ON #############')
    
    def halogen_off(self):
        self.gpio.gpio_off(7)
        if self.verbose:
            print('#### DH_mini : halogen lamp OFF #############')


    def shutter_open(self):
        self.gpio.gpio_on(6)
        if self.verbose:
            print('#### DH_mini : shutter OPEN #############')
    
    def shutter_close(self):
        self.gpio.gpio_off(6)
        if self.verbose:
            print('#### DH_mini : shutter CLOSE #############')

    def check_status(self):
        status = {
            'halogen' : self.num_to_status(int(self.gpio.gpio_read(7))),
            'deuterium' : self.num_to_status(int(self.gpio.gpio_read(5))), 
            'shutter' : self.num_to_status(int(self.gpio.gpio_read(6)))
        }
        return status    
    
    def num_to_status(self, status_num):
        if status_num == 0:
            return 'off'
        elif status_num == 1:
            return 'on'

    
    





#serial_number = APT.list_available_devices()[0][1]
#print(APT.hardware_info(serial_number))

#print(serial_number)
#print(serial_number)
#if serial_number:
#Rotor = ThorlabsK10CR1(55119384, sethome=True)
#Rotor.move_by(15)
#print(Rotor.get_position())
# Rotor._move_home()
# print(Rotor._get_position())
#Rotor.move_to(10)
#print(Rotor._get_position())
