from datetime import datetime, date
import os

import socket


HEADER_LENGTH=10

def new_log(Rxn_Name):
    Date_Tag = str(date.today())
    New_Name = os.getcwd() + '\\log_file\\' + Rxn_Name + '_' + Date_Tag + '.log'
    with open(New_Name,'w') as New_Log:
        Time_Tag = str(datetime.now())
        First_Line = Time_Tag + '\t' + Rxn_Name + ' Start Run\n'
        New_Log.write(First_Line)
    return New_Name

def write_log(File,Event):
    #get starting time in ms
    
    #open log file as append
    with open(File,'a') as Log:
        Time_Tag = str(datetime.now())
        To_Write = Time_Tag + '\t' + Event + '\n'
        print(To_Write)
        Log.write(To_Write)


class Logger():
    def __init__(
        self,
        socket_client,
    ) -> None:

        self.silasocket = socket_client
        return

    def new_log(self, Rxn_Name):
        Date_Tag = str(date.today())
        New_Name = os.getcwd() + '\\log_file\\' + Rxn_Name + '_' + Date_Tag + '.log'
        with open(New_Name,'w') as New_Log:
            Time_Tag = str(datetime.now())
            First_Line = Time_Tag + '\t' + Rxn_Name + ' Start Run\n'
            New_Log.write(First_Line)
        return New_Name

    def write_log(self, File,Event):
        #get starting time in ms
        
        #open log file as append
        with open(File,'a') as Log:
            Time_Tag = str(datetime.now())
            To_Write = Time_Tag + '\t' + Event + '\n'
            print(To_Write)
            Log.write(To_Write)
        
        message = Event.encode('utf-8')
        message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
        self.silasocket.send(message_header+ message)




    

