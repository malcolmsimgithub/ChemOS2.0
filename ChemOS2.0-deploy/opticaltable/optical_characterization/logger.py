import socket
from datetime import datetime
HEADER_LENGTH = 10
IP = "127.0.0.1"
PORT = 65012



class Logger(object):
    def __init__(self, stdout=True, logfile=None, pause = False, callback=None, time_format='({:%Y/%m/%d %H:%M:%S})'):
        self.stdout = stdout
        self.logfile = logfile
        self.pause = pause
        self.time_format = time_format
        self.callback = callback
        self.socketid = 'optics-table'.encode('utf-8')
        
        self.silasocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # Connect to a given ip and port
        self.silasocket.connect((IP, PORT))
        # Set connection to non-blocking state, so .recv() call won;t block, just return some exception we'll handle
        self.silasocket.setblocking(False)

        id  = 'optics-table'.encode('utf-8')
        username_header = f"{len(id):<{HEADER_LENGTH}}".encode('utf-8')
        self.silasocket.send(username_header + id)

    def print(self, message):
        "In message, use {timestamp} to format"
        message = message.format(timestamp=self.time_format.format(datetime.now()))

        if self.callback:
            self.callback(message)

        if self.logfile:
            with open(self.logfile, 'a') as f:
                f.write(message + '\n')
        if self.stdout:
            if self.pause:
                _ = input(f'{message} (pause)')
            else:
                print(message)
        
        message = (message).encode('utf-8')
        message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
        self.silasocket.send(message_header+ message)


    def __call__(self, message):
        self.print(message)

