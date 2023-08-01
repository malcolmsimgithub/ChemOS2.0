import drivers.valve as valve
import drivers.pump as pump
import drivers.hotplate as hotplate
#import drivers.solenoid as solenoid
import drivers.eventlog as eventlog
from drivers.eventlog import Logger
import numpy as np
import time
import sys
import os
import socket

SOCKET_IP = "127.0.0.1"
SOCKET_PORT = 65001
SOCKET_ID = 'The machine'
HEADER_LENGTH=10

N2_Port = 8
Liquid_G1 = 7
Liquid_G2 = 6
Waste_Port = 5
Shared_Valve = 4

Gas_Speed = 50
Gas_Speed_H = 150
Org_Speed = 15
Aqu_Speed = 10

RPM = 500

Module = 2

#create logger and cnnect to socket
client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect((SOCKET_IP, SOCKET_PORT))
client_socket.setblocking(False)
username = SOCKET_ID.encode('utf-8')
username_header = f"{len(username):<{HEADER_LENGTH}}".encode('utf-8')
client_socket.send(username_header + username)

logger = Logger(client_socket)

#change reaction vial atmosphere constant No.
def rxn_atm_change_const(Sample_No,Cycle_No):
    valve.valve_move(Module,Shared_Valve,8)
    for i in range (0,Sample_No):
        Rxn_Group = i // 8 + 1        #get Valve_No of sample to add to
        Rxn_Port = i % 8 + 1          #get Port_No of sample to add to
        valve.valve_move(Module,Rxn_Group,Rxn_Port)
        for j in range (0,Cycle_No):
            pump.pump_move(Module,Rxn_Group,Gas_Speed,10)         #draw 10 ml gas from vial
            pump.pump_move(Module,Waste_Port,Gas_Speed_H,0)           #exhause to waste
            print('Atmosphere changed ', j)
        Event = 'Sample ' + str(i+1) + ' Atmosphere Changed'
        logger.write_log(Log_File,Event)

#get liquid from Shared_Valve 
def rxn_get_liquid(Module,Liquid_No,Amount,Rinse=False):
    pump.pump_move(Module,Waste_Port,Org_Speed,0)
    if Liquid_No == -1:
        if Rinse:
            for i in range (0,2):
                pump.pump_move(Module,Liquid_G1,Org_Speed,0.8)
                pump.pump_move(Module,Waste_Port,Org_Speed,0)
            Event = 'Syringe ' + str(Module) + 'rinsed with Liquid G_1'
            logger.write_log(Log_File,Event)
        pump.pump_move(Module,Liquid_G1,Org_Speed,Amount)
        Event = str(Amount) + 'mL Liquid G_1 in Syringe ' + str(Module)
        logger.write_log(Log_File,Event)
    elif Liquid_No == -2:
        if Rinse:
            for i in range (0,2):
                pump.pump_move(Module,Liquid_G2,Org_Speed,0.8)
                pump.pump_move(Module,Waste_Port,Org_Speed,0)
            Event = 'Syringe ' + str(Module) + 'rinsed with Liquid G_1'
            logger.write_log(Log_File,Event)
        pump.pump_move(Module,Liquid_G2,Org_Speed,Amount)
        Event = str(Amount) + 'mL Liquid G_1 in Syringe ' + str(Module)
        eventlog.write_log(Log_File,Event)
    elif 31 <= Liquid_No <= 38:      
        valve.valve_move(Module,3,(Liquid_No - 30))        
        if Rinse:
            for i in range (0,3):
                pump.pump_move(Module,3 ,Org_Speed,0.5)
                pump.pump_move(Module,Waste_Port,Org_Speed,0)
            Event = 'Syringe ' + str(Module) + 'rinsed with liquid ' + str(Liquid_No)
            logger.write_log(Log_File,Event)
        pump.pump_move(Module,3 ,Org_Speed,Amount)
        Event = str(Amount) + 'mL Liquid ' + str(Liquid_No) + ' in Syringe ' + str(Module)
        logger.write_log(Log_File,Event)
    elif 41 <= Liquid_No <= 48:      
        valve.valve_move(Module,4,(Liquid_No - 40))        
        if Rinse:
            for i in range (0,3):
                pump.pump_move(Module,4,Org_Speed,0.5)
                pump.pump_move(Module,Waste_Port,Org_Speed,0)
            Event = 'Syringe ' + str(Module) + 'rinsed with liquid ' + str(Liquid_No)
            logger.write_log(Log_File,Event)
        pump.pump_move(Module,4,Org_Speed,Amount)
        Event = str(Amount) + 'mL Liquid ' + str(Liquid_No) + ' in Syringe ' + str(Module)
        logger.write_log(Log_File,Event)

#add liquid to Rxn
def rxn_add_liquid(Module,Rxn_Group,Rxn_Port,Flush=True):
    valve.valve_move(Module,Rxn_Group,Rxn_Port)
    pump.pump_move(Module,Rxn_Group,Org_Speed,0)
    Event = 'Liquid add to Group ' + str(Rxn_Group) + ' Rxn_Port ' + str(Rxn_Port)
    logger.write_log(Log_File,Event)
    if Flush:
        for i in range (0,3):
            pump.pump_move(Module,N2_Port,Gas_Speed,2)
            pump.pump_move(Module,Rxn_Group,Gas_Speed,0)
        Event = 'Tubing flushed with N2'
        logger.write_log(Log_File,Event)

#make stock solution
def make_stock_solution(Module, Stock_Group, Stock_Port, Volume):
    rxn_get_liquid(Module,Liquid_No = -1, Amount = 0,Rinse = True)  #pre-rinse syringe
    Liq = Stock_Group * 10 + Stock_Port
    valve.valve_move(Module,Stock_Group,Stock_Port)
    print('Change atmosphere of stock on ', Liq)
    for i in range (0,10):                                          #change atmosphere of stock solution vial
        pump.pump_move(Module,Stock_Group,Gas_Speed,10)             #draw 10 ml gas from vial
        pump.pump_move(Module,Waste_Port,Gas_Speed_H,0)               #exhause to waste
        print('Atmosphere changed ', i)
    Event = 'Stock ' + str(Stock_Group) + str(Stock_Port) + ' Atmosphere Changed'
    logger.write_log(Log_File,Event)
    rxn_get_liquid(Module,Liquid_No = -1, Amount = Volume,Rinse = True) #add pure solvent to stock solution vial
    rxn_add_liquid(Module, Stock_Group, Stock_Port, Flush = True)
    for j in range (0, 5):                                              #mix well
        rxn_get_liquid(Module, Liquid_No = Liq, Amount = Volume, Rinse = False)
        rxn_add_liquid(Module, Stock_Group, Stock_Port, Flush = True)
        print('Solution mixed ', j)
    rxn_get_liquid(Module, Liquid_No = Liq, Amount = Volume, Rinse = False)
    rxn_add_liquid(Module, Stock_Group, Stock_Port, Flush = False)
    Event = 'Stock ' + str(Stock_Group) + str(Stock_Port) + ' Preparation finished'
    logger.write_log(Log_File,Event)
    
    
'''
#flush with N2 flow
def rxn_N2_flush(Solenoid_No,Time=600):
    solenoid.solen_operation(Solenoid_No,1)
    Event = 'N2 flush starts'
    eventlog.write_log(Log_File,Event)
    time.sleep(Time)
    solenoid.solen_operation(Solenoid_No,0)
    Event = 'Vials ' + str(Solenoid_No + 1) + 'flushed with N2 for ' + str(Time) + 's'
    eventlog.write_log(Log_File,Event)
'''     
#heat and stir
def rxn_heat_stir(Temp,RPM,Time):
    hotplate.hotplate_stir(True,RPM)
    hotplate.hotplate_temp(True,Temp)
    Event = 'Start Heating ' + str(Temp) + ' Stirring ' + str(RPM) + ' for ' + str(Time) + ' hours'
    logger.write_log(Log_File,Event)                                          #write event to log file
    time.sleep(Time*3600)
    hotplate.hotplate_stir(False)
    hotplate.hotplate_temp(False)
    Event = 'Stop Heating and Stirring'
    logger.write_log(Log_File,Event)

#input reaction name, and find if setting parameter .csv file exist

if len(sys.argv)== 2:
    Rxn_Name = sys.argv[1]
else:
    print('Please input reaction name to save log file as:')
    Rxn_Name = input()
Setting_File = Rxn_Name + '.csv'
with open(Setting_File,'r',encoding='utf-8-sig') as Rxn_Setting:
    Log_File = eventlog.new_log(Rxn_Name)


    
#read settings, Sample_No, Cycle_No, Temperature, Rxn Time
with open(Setting_File,'r',encoding='utf-8-sig') as Rxn_Setting:
    for i, line in enumerate(Rxn_Setting):
        if i == 1:
            line = line.rstrip('\n')
            Setting_List = line.split(',')
            Sample_No = int(Setting_List[0])
            Cycle_No = int(Setting_List[1])
            Flush_Time = int(Setting_List[2])
            #RPM = int(Setting_List[3])
        elif i == 3:
            line = line.rstrip('\n')
            Setting_List = line.split(',')
            Rxn_Temp_1 = int(Setting_List[0])
            Rxn_Time_1 = int(Setting_List[1])
            Rxn_Temp_2 = int(Setting_List[2])
            Rxn_Time_2 = int(Setting_List[3])
        elif i > 3:
            break

        
#initialize pump
if pump.pump_initialize():
    Event = 'Pump Initialization Successful'
    logger.write_log(Log_File,Event)                                            #write event to log file
else:
    Event = 'Pump Initialization Failed'
    logger.write_log(Log_File,Event)                                            #write event to log file
    sys.exit()

#read liquid addition parameters, in format of Liquid_No, Destination, Flush tubing (T/F)
with open(Setting_File,'r',encoding='utf-8-sig') as Liquid_Amount:
    Liquid_Amount_Array = np.genfromtxt(Liquid_Amount,dtype=float,delimiter=',',skip_header=5)

#change reaction atmosphere
if Cycle_No:
    rxn_atm_change_const(Sample_No,Cycle_No)
#if Flush_Time:
#    rxn_N2_flush(Solenoid_No=3,Time=(Flush_Time*60))


#add liquids to reaction
#Total_Addition = np.shape(Liquid_Amount_Array)[0]
for i in range (0,16):
    Rxn_Group = int(Liquid_Amount_Array[i,1]) // 10         #get Valve_No of sample to add to
    Rxn_Port = int(Liquid_Amount_Array[i,1]) % 10          #get Port_No of sample to add to
    Flush = bool(Liquid_Amount_Array[i,3])                      #get Flush tubing after addition (True or False)
    if ((i == 0) or (Liquid_Amount_Array[i,0] != Liquid_Amount_Array[(i-1),0])):                #rinse tubing and syring for first sample of every solution
        rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),0,Rinse=True)
    if Liquid_Amount_Array[i,2]:
        rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),Liquid_Amount_Array[i,2],Rinse=False)
        rxn_add_liquid(Module,Rxn_Group,Rxn_Port,Flush)


#start heating and stirring
if Rxn_Time_1:
    rxn_heat_stir(Rxn_Temp_1,RPM,Rxn_Time_1)
   
#add liquids to reaction
Total_Addition = np.shape(Liquid_Amount_Array)[0]
for i in range (16,32):
    Stock_Group = int(Liquid_Amount_Array[i,0]) // 10         #get Valve_No of sample to add to
    Stock_Port = int(Liquid_Amount_Array[i,0]) % 10          #get Port_No of sample to add to
    
    Rxn_Group = int(Liquid_Amount_Array[i,1]) // 10         #get Valve_No of sample to add to
    Rxn_Port = int(Liquid_Amount_Array[i,1]) % 10          #get Port_No of sample to add to
    Flush = bool(Liquid_Amount_Array[i,3])                      #get Flush tubing after addition (True or False)
    #print(Stock_Group, Stock_Port)
    if ((i == 0) or (Liquid_Amount_Array[i,0] != Liquid_Amount_Array[(i-1),0])):                #rinse tubing and syring for first sample of every solution
        make_stock_solution(Module, Stock_Group, Stock_Port, Volume = 6)
        #rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),0,Rinse=False)
    if Liquid_Amount_Array[i,2]:
        rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),Liquid_Amount_Array[i,2],Rinse=False)
        rxn_add_liquid(Module,Rxn_Group,Rxn_Port,Flush)

for i in range (32,Total_Addition):
    Rxn_Group = int(Liquid_Amount_Array[i,1]) // 10         #get Valve_No of sample to add to
    Rxn_Port = int(Liquid_Amount_Array[i,1]) % 10          #get Port_No of sample to add to
    Flush = bool(Liquid_Amount_Array[i,3])                      #get Flush tubing after addition (True or False)
    if ((i == 0) or (Liquid_Amount_Array[i,0] != Liquid_Amount_Array[(i-1),0])):                #rinse tubing and syring for first sample of every solution
        #make_stock_solution(Module, Rxn_Group, Rxn_Port, Volume = 6)
        rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),0,Rinse=True)
    if Liquid_Amount_Array[i,2]:
        rxn_get_liquid(Module,int(Liquid_Amount_Array[i,0]),Liquid_Amount_Array[i,2],Rinse=False)
        rxn_add_liquid(Module,Rxn_Group,Rxn_Port,Flush)
        
#start heating and stirring
if Rxn_Time_2:
    rxn_heat_stir(Rxn_Temp_2,RPM,Rxn_Time_2)

hotplate.hotplate_close()
#solenoid.solen_close()


