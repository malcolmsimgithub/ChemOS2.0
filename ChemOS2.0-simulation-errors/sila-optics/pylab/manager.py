from . import instruments as ins
import numpy as np
import matplotlib.pyplot as plt
import time, json
from datetime import datetime
import socket

HEADER_LENGTH = 10
IP = "127.0.0.1"
PORT = 65012


class Manager(object):
    def __init__(self, settings):
        if isinstance(settings, str):
            with open(settings, 'r') as f:
                settings = json.load(f)
        elif not isinstance(settings, list):
            raise ValueError('Settings should be filename or list')

        self.devices = {}

        # register equipments
        for instrument_settings in settings:
            self.register_instrument(instrument_settings)

    def register_instrument(self, settings):
        i_name = settings.get('name', None)
        i_type = settings.get('instrument', None)
        i_sett = settings.get('settings', {})

        if i_name is None or i_type is None:
            raise ValueError('Each instrument settings should have name and instrument type (instrument)')

        # register instrument
        print(f'Register instrument {i_type} as {i_name}')
        self.devices[i_name] = getattr(ins, i_type)(**i_sett)

    def __getattr__(self, attr):
        device = self.devices.get(attr, None)
        if device:
            return device
        else:
            raise ValueError('No instrument')



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

