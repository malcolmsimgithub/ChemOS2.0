import socket

HEADER_LENGTH = 10
IP = "127.0.0.1"
PORT = 65005

class Logger():
    def __init__(self) -> None:

        self.client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Connect to a given ip and port
        self.client_socket.connect((IP, PORT))

        # Set connection to non-blocking state, so .recv() call won;t block, just return some exception we'll handle
        self.client_socket.setblocking(False)

        id  = 'HPLCMS'.encode('utf-8')
        username_header = f"{len(id):<{HEADER_LENGTH}}".encode('utf-8')
        self.client_socket.send(username_header + id)

    def log(self, message):

        print(message)
        # Prepare username and header and send them
        # We need to encode username to bytes, then count number of bytes and prepare header of fixed size, that we encode to bytes as well
        message = message.encode('utf-8')
        message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
        self.client_socket.send(message_header + message)
