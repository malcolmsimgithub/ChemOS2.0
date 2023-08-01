import visa
import time

#communicate to devices and list them
rm = visa.ResourceManager()

#connect to desired COM port device and open remote control
#COM_No = 6
COM_No = 4
COM_Port = 'ASRL'+str(COM_No)+'::INSTR'
hotplate = rm.open_resource(COM_Port)
simulation = True

#set temperature 
def hotplate_temp(Heat_ONOFF=False,Temp=20):
    if simulation:
        time.sleep(1)
        return None

    if Heat_ONOFF == True:
        #set temp
        hotplate.write('OUT_SP_1 '+str(Temp))
        time.sleep(1)
        #start heating
        hotplate.write('START_1')
        time.sleep(1)
    else:
        #stop heating
        hotplate.write('STOP_1')
        time.sleep(1)

#set stirring
def hotplate_stir(Stir_ONOFF=False,RPM=0):
    if simulation:
        time.sleep(1)
        return None

    if Stir_ONOFF == True:
        #set rpm
        hotplate.write('OUT_SP_4 '+str(RPM))
        time.sleep(1)
        #start stirring
        hotplate.write('START_4')
        time.sleep(1)
    else:
        #stop stirring
        hotplate.write('STOP_4')
        time.sleep(1)

#weigh mass
def hotplate_weigh(Tare_ONOFF):
    if simulation:
        return None

    if Tare_ONOFF == True:
        #reset taring value
        hotplate.write('STOP_90')
        hotplate.write('START_90')
        time.sleep(10)
        hotplate.write('STATUS_90')
        print(hotplate.read())
        return 0
    else:
        #check stability
        for i in range (0,6):
            hotplate.write('STATUS_90')
            Hotplate_Reading = hotplate.read()
            Hotplate_Reading = Hotplate_Reading.strip()
            if Hotplate_Reading=='1041 90':
            #    print('y')
                break
            else:
             #   print('n')
                time.sleep(10)
            #out put error here?
        #measure weight
        hotplate.write('IN_PV_90')
        time.sleep(1)
        Weight = hotplate.read()
        return Weight

#close control, do NOT shut down hotplate!!!
def hotplate_close():
    hotplate.close()

#hotplate_stir(False,300)
#hotplate_temp(False)
