import visa
import time
#COM_List=[10,7,11,10,9,9,8]
COM_List=[9,9,9,9,9,9]
rm = visa.ResourceManager()

def valve_move(Valve_Module,Valve_No,Valve_Port):
    #select COM Port to use
    COM_No = COM_List[Valve_Module]
    #print(COM_No)
    COM_Port = 'ASRL'+str(COM_No)+'::INSTR'
    valve = rm.open_resource(COM_Port)
    #build command, then move valve
    Command = '/'+str(Valve_No)+'o'+str(Valve_Port)+'R'
    valve.write(Command)
    time.sleep(2)
    valve.write(Command)
    time.sleep(2)
    valve.close()
    #time.sleep(2)

'''
time.sleep(5)
valve_move(1,1,8)
time.sleep(5)
valve_move(1,1,7)
time.sleep(5)
valve_move(1,1,5)
time.sleep(5)
valve_move(1,1,4)
'''

#valve_move(5, 1, 5)
