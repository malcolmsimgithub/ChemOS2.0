import socket
import os
from sila2.server import MetadataDict, ObservableCommandInstanceWithIntermediateResponses
import socket

from ..generated.hplcmssimulator import (
    GetResults1st_IntermediateResponses,
    GetResults1st_Responses,
    HPLCMSsimulatorBase,
    Status_Responses,
    SubmitJobAutosampler_IntermediateResponses,
    SubmitJobAutosampler_Responses,
    SubmitJobChemspeed_IntermediateResponses,
    SubmitJobChemspeed_Responses,
    ValveStatus_Responses,
)

from pathlib import Path
import sys
import errno
import pickle
import subprocess


## coordinates for the socket server connection
IP = "127.0.0.1"
PORT = 65005
SOCKET_ID = 'sila2'
HEADER_LENGTH = 10

filepath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

JOBS_RETURNED = os.path.abspath("data_path/jobs_returned/")


def get_status(statusfile):
    with open(statusfile, "rb") as f:
        status = pickle.load(f)['status']
    return status

def run_chemspeed_client(client_socket, instance, outputfile, jobname):

    outputfile = os.path.join(JOBS_RETURNED, f"{jobname}.json")
    outputdir = os.path.join(JOBS_RETURNED, jobname)


    bindata = str.encode("dummy")

    while True:
        try:
            while True:

                 # check for output data
                if os.path.isfile(outputfile):

                    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("job termination detected", bindata))

                    outputfile7z = outputfile[:-5]+".7z"

                    save_cmd = f"7z a {outputfile7z} {outputdir}"

                    compressjob = subprocess.run(save_cmd,shell=True,stdout=subprocess.PIPE)

                    if compressjob.returncode == 0:
                        with open(f"{outputfile7z}", "rb") as f:
                            data = f.read()
                        
                        instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("output_data", data))
                    
                    else:
                        instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("send return data failed", data))
                    return   "Characterization complete"

                username_header = client_socket.recv(HEADER_LENGTH)

                if not len(username_header):
                    print('Connection closed by the server')
                    sys.exit()

                # Convert header to int value
                username_length = int(username_header.decode('utf-8').strip())

                # Receive and decode username
                username = client_socket.recv(username_length).decode('utf-8')

                # Now do the same for message (as we received username, we received whole message, there's no need to check if it has any length)
                message_header = client_socket.recv(HEADER_LENGTH)
                message_length = int(message_header.decode('utf-8').strip())
                message = client_socket.recv(message_length).decode('utf-8')

                # Print message
                print(f'{username} > {message}')
                instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses(message, bindata))
                if message == f"job {jobname} completed" or os.path.isfile(outputfile):
                    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("job termination detected", bindata))

                    outputfile7z = outputfile[:-5]+".7z"

                    save_cmd = f"7z a {outputfile7z} {outputdir}"

                    compressjob = subprocess.run(save_cmd,shell=True,stdout=subprocess.PIPE)

                    if compressjob.returncode == 0:
                        with open(f"{outputfile7z}", "rb") as f:
                            data = f.read()
                        
                        instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("output_data", data))

                        time.sleep(30)
                    
                    else:
                        instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("send return data failed", data))
                    return   "Characterization complete"
        



        except IOError as e:
            # This is normal on non blocking connections - when there are no incoming data error is going to be raised
            # Some operating systems will indicate that using AGAIN, and some using WOULDBLOCK error code
            # We are going to check for both - if one of them - that's expected, means no incoming data, continue as normal
            # If we got different error code - something happened
            if e.errno != errno.EAGAIN and e.errno != errno.EWOULDBLOCK:
                print('Reading error: {}'.format(str(e)))
                sys.exit()

            # We just did not receive anything
            continue

        except Exception as e:
            # Any other exception - something happened, exit
            print('Reading error: '.format(str(e)))
            sys.exit()

def run_autosampler_client(client_socket, instance, outputfile, jobname):
    bindata = str.encode("dummy")

    outputfile = os.path.join(JOBS_RETURNED, f"{jobname}.json")
    outputdir = os.path.join(JOBS_RETURNED, jobname)

    
    instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses(f"monitoring for {outputfile}", bindata))

    while True:
        try:
            while True:

                 # check for output data
                if os.path.isfile(outputfile):

                    ###### TODO: convert output file to compressed 7z format to send to ifog

                    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("job termination detected", bindata))

                    outputfile7z = outputfile[:-5]+".7z"

                    save_cmd = f"7z a {outputfile7z} {outputdir}"

                    compressjob = subprocess.run(save_cmd,shell=True,stdout=subprocess.PIPE)

                    if compressjob.returncode == 0:
                        with open(f"{outputfile7z}", "rb") as f:
                            data = f.read()
                        
                        instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses("output_data", data))

                        time.sleep(30)
                    
                    else:
                        instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses("send return data failed", data))
                    return   "Characterization complete"

                username_header = client_socket.recv(HEADER_LENGTH)

                if not len(username_header):
                    print('SiLA2: Connection closed by the server')
                    sys.exit()

                # Convert header to int value
                username_length = int(username_header.decode('utf-8').strip())

                # Receive and decode username
                username = client_socket.recv(username_length).decode('utf-8')

                # Now do the same for message (as we received username, we received whole message, there's no need to check if it has any length)
                message_header = client_socket.recv(HEADER_LENGTH)
                message_length = int(message_header.decode('utf-8').strip())
                message = client_socket.recv(message_length).decode('utf-8')

                # Print message
                print(f'{username} > {message}')
                instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses(message, bindata))
                if message == f"job {jobname} completed" or os.path.isfile(outputfile):
                    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("job termination detected", bindata))

                    outputfile7z = outputfile[:-5]+".7z"

                    save_cmd = f"7z a {outputfile7z} {outputdir}"

                    compressjob = subprocess.run(save_cmd,shell=True,stdout=subprocess.PIPE)

                    if compressjob.returncode == 0:
                        with open(f"{outputfile7z}", "rb") as f:
                            data = f.read()
                        
                        instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses("output_data", data))

                        time.sleep(30)
                    
                    else:
                        instance.send_intermediate_response(SubmitJobAutosampler_IntermediateResponses("send return data failed", data))
                    return   "Characterization complete"
        

        except IOError as e:
            # This is normal on non blocking connections - when there are no incoming data error is going to be raised
            # Some operating systems will indicate that using AGAIN, and some using WOULDBLOCK error code
            # We are going to check for both - if one of them - that's expected, means no incoming data, continue as normal
            # If we got different error code - something happened
            if e.errno != errno.EAGAIN and e.errno != errno.EWOULDBLOCK:
                print('Reading error: {}'.format(str(e)))
                sys.exit()

            # We just did not receive anything
            continue

        except Exception as e:
            # Any other exception - something happened, exit
            print('Reading error: '.format(str(e)))
            sys.exit()

def update_valve_status(status):
        """
        Update the status file of the 6-port valve and the flow selector. Filename is specified in the setting. 

        Parameters
        --------------------
        Status: str
            New status of the valve ('free' or 'busy')
        """  
        fname = Path('data_path/status/valve_status_HPLCMS.pkl')
        with open (fname, 'wb') as f:
            pickle.dump(status, f)

        
        return "succesfully updated valve status"


def check_valve_status():
    """
    Check the status of the 6-port valve and the flow selector by reading the file specified in the setting. 

    Returns
    --------------------
    status: str
        Status of the valve ('free' or 'busy')
    """ 
    fname = ('data_path/status/valve_status_characterization.pkl')
    with open (fname, 'rb') as f:
        status = pickle.load(f)
        print('characterization_valve_status : %s' %status)
    return status


import time


def timestamp_date():
    """Get a timestamp string in the YYYY-MM-DD format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%y-%m-%d", time.localtime())
    return timestamp


def timestamp_time():
    """Get a timestamp string in the HH-MM format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%H-%M", time.localtime())
    return timestamp


def timestamp_time_precise():
    """Get a timestamp string in the HH-MM-SS format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime("%H-%M-%S", time.localtime())
    return timestamp


def timestamp_datetime():
    """Get a timestamp string in the YYYY-MM-DD_HH-MM format.

    Parameters:

    Returns:
        timestamp (str)
    """
    timestamp = time.strftime(f"{timestamp_date()}_{timestamp_time_precise()}", time.localtime())
    return timestamp
