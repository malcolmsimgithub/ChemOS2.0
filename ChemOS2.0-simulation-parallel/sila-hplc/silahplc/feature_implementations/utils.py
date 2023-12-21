import os
from ..generated.hplcmssimulator import (
    SubmitJobChemspeed_IntermediateResponses,
    SubmitJobChemspeed_Responses,
    BlankRun_IntermediateResponses,
    LostConnection
)


from sila2.framework import SilaError

from pathlib import Path
import sys
import errno
import pickle
import json
import time

## coordinates for the socket server connection
IP = "127.0.0.1"
PORT = 65005
SOCKET_ID = 'sila2'
HEADER_LENGTH = 10

filepath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

JOBS_RETURNED = os.path.abspath("communication/jobs_returned/")



def get_status(statusfile):
    with open(statusfile, "rb") as f:
        status = pickle.load(f)['status']
    return status

def run_chemspeed_client(client_socket, instance, outputfile, jobname):

    outputfile = os.path.join(JOBS_RETURNED, f"{jobname}.json")
    bindata = str.encode("dummy")
    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses(f"monitoring for {outputfile}", bindata))

    outputfile7z = outputfile[:-5]+".7z"
    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses(f"monitoring for {outputfile7z}", bindata))

    while True:
        try:
            while True:
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
                    outputjson = outputfile[:-5]+".json"
                    time.sleep(10)

                    with open(f"{outputfile7z}", "rb") as f:
                        data = f.read()
                    
                    instance.send_intermediate_response(SubmitJobChemspeed_IntermediateResponses("output_data", data))

                    time.sleep(30)

                    with open(outputjson, "r") as f:
                        output_json = json.load(f)
                    
                    return  str(output_json["result"])
                
                elif message == "HPLCMS_lost":
                    raise LostConnection

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

def run_blank_client(client_socket, instance, outputfile, jobname):

    outputfile = os.path.join(JOBS_RETURNED, f"{jobname}.json")

    bindata = str.encode("dummy")
    outputfile7z = outputfile[:-5]+".7z"
    instance.send_intermediate_response(BlankRun_IntermediateResponses(f"monitoring for {outputfile7z}", bindata))

    while True:
        try:
            while True:
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
                instance.send_intermediate_response(BlankRun_IntermediateResponses(message, bindata))
                if message == f"job {jobname} completed" or os.path.isfile(outputfile):
                    instance.send_intermediate_response(BlankRun_IntermediateResponses("job termination detected", bindata))

                    outputfile7z = outputfile[:-5]+".7z"

                    time.sleep(10)

                    with open(f"{outputfile7z}", "rb") as f:
                        data = f.read()
                    
                    instance.send_intermediate_response(BlankRun_IntermediateResponses("output_data", data))

                    time.sleep(30)
                
                    return   "Run complete"
            

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


