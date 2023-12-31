# Generated by sila2.code_generator; sila2.__version__: 0.10.1
from __future__ import annotations

from typing import TYPE_CHECKING

from sila2.server import MetadataDict, ObservableCommandInstanceWithIntermediateResponses

from ..generated.chemspeedoperator import (
    Addbatch_IntermediateResponses,
    Addbatch_Responses,
    AddCharacterization_IntermediateResponses,
    AddCharacterization_Responses,
    ChangeProcedure_Responses,
    ChemSpeedOperatorBase,
)

import os
import shutil
import subprocess
from chemspeed_operator_process.Utils import timestamp_datetime
import time
from typing import TYPE_CHECKING
import socket
import sys
import errno
from pathlib import Path
import pandas as pd
import json
import glob
if TYPE_CHECKING:
    from ..server import Server


from .utils import isbusy, run_synthesis_client, run_characterization_client

OUTPUT: Path = Path("chmspd_output")
STORAGE: Path =Path('chmspd_storage')
BATCHFOLDER: Path = os.path.join(os.path.join(OUTPUT, "Synthesis"), "Batches_to_Make")
PROCEDUREDIR = Path(os.path.join(OUTPUT, "Synthesis"))/ "Procedure"
COMPLETEDDIR = os.path.join(os.path.join(OUTPUT, "Synthesis"), "Completed_Batches")
CHARFILE = os.path.join(os.path.join(OUTPUT, "Characterization"), "Positions_for_Characterization.csv")
HEADER_LENGTH = 10
IP = "127.0.0.1"
PORT = 65041
SOCKET_ID = 'sila2'
my_username = "sila2_synthesis"

class ChemSpeedOperatorImpl(ChemSpeedOperatorBase):
    def __init__(self, parent_server: Server) -> None:
        super().__init__(parent_server=parent_server)

    def get_IsBusy(self, *, metadata: MetadataDict) -> str:

        

        return isbusy()

    def ChangeProcedure(self, SynthesisProcedure: str, *, metadata: MetadataDict) -> ChangeProcedure_Responses:
        # if isbusy() != True:
            for file in os.listdir(PROCEDUREDIR):
                file_path = PROCEDUREDIR/file
                os.remove(file_path)
            shutil.copyfile(STORAGE/SynthesisProcedure, PROCEDUREDIR/SynthesisProcedure)
            return ChangeProcedure_Responses("procedure succesfully changed")
        # else:
        #     return ChangeProcedure_Responses("Chemspeed is currently busy")
        

    def Addbatch(
        self,
        BatchName: str,
        Batchfile: str,
        *,
        metadata: MetadataDict,
        instance: ObservableCommandInstanceWithIntermediateResponses[Addbatch_IntermediateResponses],
    ) -> Addbatch_Responses:
        instance.begin_execution() 

        evaldict = eval(Batchfile) 

        datetime = timestamp_datetime()
        batchfile = f"{BatchName}_{datetime}.json"

        batchtoken = batchfile.strip(".json")
        with open( os.path.join(BATCHFOLDER,batchfile), 'w') as f:
            json.dump(evaldict, f)
        instance.send_intermediate_response(Addbatch_IntermediateResponses("starting synthesis", "operation string"))


        # Create a socket
        client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        client_socket.connect((IP, PORT))
        client_socket.setblocking(False)
        username = my_username.encode('utf-8')
        username_header = f"{len(username):<{HEADER_LENGTH}}".encode('utf-8')
        client_socket.send(username_header + username)
        return Addbatch_Responses(run_synthesis_client(client_socket, instance, batchtoken, batchfile))
         
    def AddCharacterization(
        self,
        RackPosition: str,
        Identifier: str,
        Structure: str,
        FilterCollect: str,
        *,
        metadata: MetadataDict,
        instance: ObservableCommandInstanceWithIntermediateResponses[AddCharacterization_IntermediateResponses],
    ) -> AddCharacterization_Responses:
        instance.begin_execution()  # set execution status from `waiting` to `running`

        df = pd.DataFrame

        jobname = f"{Identifier}_{timestamp_datetime()}"
    
        jobsfile = "/Users/maozer/VSCodeProjects/sila-chemspeed/chmspd/output/Characterization/Characterizations.json"

        datadict = {}

        datadict['Identifier'] = jobname
        datadict['Structure'] = Structure
        datadict['filter_collect'] = bool(FilterCollect)
        datadict['characterization_1st']= True
        datadict['characterization_2nd']= False

        HEADER_LENGTH = 10
        IP = "127.0.0.1"
        PORT = 65001
        my_username = "sila2_characterization"

        # Create a socket
        client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        client_socket.connect((IP, PORT))
        client_socket.setblocking(False)
        username = my_username.encode('utf-8')
        username_header = f"{len(username):<{HEADER_LENGTH}}".encode('utf-8')
        client_socket.send(username_header + username)
        submission_file = pd.read_csv(CHARFILE, index_col=0)

        for column in submission_file.columns.values:
            submission_file.at[RackPosition, column] = datadict[column]
        submission_file.to_csv(CHARFILE)

        instance.send_intermediate_response(AddCharacterization_IntermediateResponses(jobname, "operations"))

        return AddCharacterization_Responses(run_characterization_client(client_socket, instance, jobname, jobsfile))  

