from thorlabs_K10CR1 import Rotor
import thorlabs_APT as APT

serial_number = APT.list_available_devices()[0][1]

if serial_number:
    print(serial_number)
    Rotor = Rotor(serial_number)
    Rotor._move_by(50)
    Rotor._move_home()
    Rotor._move_to(0)
