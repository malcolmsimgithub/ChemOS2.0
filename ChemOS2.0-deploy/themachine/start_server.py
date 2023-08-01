import os
import subprocess
from pathlib import Path
import threading
import os
from run_socket_server import run_socket_server




#################################
# INITIALIZE SILA SERVER        #
#################################

server_cmd = f"python -m silathemachine --ip-address 127.0.0.1 -p 65000 --insecure"
proc = subprocess.Popen(server_cmd,shell=True,stdout=subprocess.PIPE)


#################################
# INITIALIZE SOCKET SERVER      #
#################################
socket_thread = threading.Thread(target=run_socket_server)
socket_thread.start()
