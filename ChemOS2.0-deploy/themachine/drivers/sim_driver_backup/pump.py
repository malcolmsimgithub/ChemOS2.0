import sys
import clr
import os
import time


clr.AddReference(r"D:\Aspuru-Guzik Lab Dropbox\Lab Manager Aspuru-Guzik\PythonScript\Han\TheMachine_Sim\TheMachine_Sim\drivers\KEMPumpDLL")

#load syringe pump control commands
from KEMPumpDLL import SyringePumpDef

pump = SyringePumpDef()
JKPump_List = [0,0,]#0,0,0,0,0]
JKPump_Active = JKPump_List
Err = ""
simulation = True
#initialize pump,out put online pump to JKPump_List as 1
def pump_initialize():
    
    #open communication
    #if pump.OpenCommunications():
    if simulation: 
        #discover module
        for i in range (0,len(JKPump_List)):
            #if pump.DiscoverModule(i+1):
            if simulation:
                #print("Pump Found No.",(i+1))
                JKPump_List[i] = 1
        #initialize pump
        for j in range (0,len(JKPump_List)):
            if JKPump_List[j] == 1:
                #if pump.Initialize((j+1),5):
                if simulation:
                    JKPump_Active[j] =1
                    #print("Pump Initialized No.",(j+1))
                else:
                    JKPump_Active[j] =0
                    print("Pump Initialize Failed No.",(j+1))
        return True
        
    else:
        print("Communication failed")
        #return False

#set pump X to access port Y and move to volumn Z, Wait_Ready=False only for moving two pump simultaneously (BE CAREFUL!)
def pump_move(Pump_No,Pump_Port,Topspeed,Volume,Wait_Ready=True):
    if simulation:
        time.sleep(Volume)
        return None
    if JKPump_Active[(Pump_No-1)] == 0:
        Err = "Pump not active"
        print(Err)
    else:
        pump.Port(Pump_No,Pump_Port,True)
        pump.Speed(Pump_No,Topspeed)
        Position = Volume*307200
        pump.MoveToPosition(Pump_No,Position,Wait_Ready)
        
def pump_port(Pump_No,Pump_Port):
    if simulation:
        time.sleep(1)
        return None
    if JKPump_Active[(Pump_No-1)] == 0:
        Err = "Pump not active"
        print(Err)
    else:
        pump.Port(Pump_No,Pump_Port,True)

def pump_position(Pump_No):
    if simulation:
        time.sleep(1)
        return None
    if JKPump_Active[(Pump_No-1)] == 0:
        Err = "Pump not active"
        print(Err)
    else:
        Position = pump.SyringePosition(Pump_No)
        return Position

def pump_set_home(Pump_No):
    if simulation:
        time.sleep(30)
        return None
    pump.SetHome_Step1(Pump_No)
    input('Press Enter to continue\n')
    pump.SetHome_Step2(Pump_No)
    
#close control
def pump_close():
    pump.CloseCommPort()


#pump_initialize()
#pump_move(2,5,20,10)
#input()
#pump_move(2,3,20,9.5)
#input()
#pump_move(2,4,20,9)
#pump_set_home(Pump_No=4)
#print(JKPump_List)
#print(JKPump_Active)
#pump_move(Pump_No = 4,Pump_Port = 5,Topspeed = 10,Volume = 2,Wait_Ready=True)
'''
for i in range (0,3):
    pump_move(Pump_No=4,Pump_Port=7,Topspeed=20,Volume=2,Wait_Ready=True)
    pump_move(Pump_No=4,Pump_Port=5,Topspeed=20,Volume=0,Wait_Ready=True)

pump_move(1,5,20,2)
pump_move(1,5,20,0)

print(pump_position(1))
pump_move(1,1,20,1)
print(pump_position(1))
pump_move(1,1,20,0)

pump_close()
'''
