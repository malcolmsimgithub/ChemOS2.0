HEADER_LENGTH = 10

IP = "127.0.0.1"
PORT = 65051


import socket
import select

def read_message(client_socket):

    try:
        # Receive our "header" containing message length, it's size is defined and constant
        message_header = client_socket.recv(HEADER_LENGTH)
        
        # If we received no data, client gracefully closed a connection
        if not len(message_header):
            return False

        # Convert header to int value
        message_length = int(message_header.decode('utf-8').strip())

        # Return an object of message header and message data
        return {'header': message_header, 'data': client_socket.recv(message_length)}

    except:
        return False

def run_socket_server():

    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server_socket.bind((IP, PORT))

    # This makes server listen to new connections
    server_socket.listen()

    # List of sockets for select.select()
    sockets_list = [server_socket]

    # List of connected clients - socket as a key, user header and name as data
    clients = {}

    print(f'Socket is listening for connections on {IP}:{PORT}...')
    while True:
        # Calls Unix select() system call or Windows select() WinSock call with three parameters:
        #   - rlist - sockets to be monitored for incoming data
        #   - wlist - sockets for data to be send to (checks if for example buffers are not full and socket is ready to send some data)
        #   - xlist - sockets to be monitored for exceptions (we want to monitor all sockets for errors, so we can use rlist)
        # Returns lists:
        #   - reading - sockets we received some data on (that way we don't have to check sockets manually)
        #   - writing - sockets ready for data to be send thru them
        #   - errors  - sockets with some exceptions
        # This is a blocking call, code execution will "wait" here and "get" notified in case any action should be taken
        read_sockets, _, exception_sockets = select.select(sockets_list, [], sockets_list)


        # Iterate over notified sockets
        for notified_socket in read_sockets:

            # If notified socket is a server socket - new connection, accept it
            if notified_socket == server_socket:

                # Accept new connection
                # That gives us new socket - client socket, connected to this given client only, it's unique for that client
                # The other returned object is ip/port set
                client_socket, client_address = server_socket.accept()

                # Client should send his name right away, receive it
                id = read_message(client_socket)

                # If False - client disconnected before he sent his name
                if id is False:
                    continue

                # Add accepted socket to select.select() list
                sockets_list.append(client_socket)

                # Also save username and username header
                clients[client_socket] = id

                print('SOCKET: Accepted new connection from {}:{}, id: {}'.format(*client_address, id['data'].decode('utf-8')))

            # Else existing socket is sending a message
            else:
                # Receive message
                message = read_message(notified_socket)

                # If False, client disconnected, cleanup
                if message is False:
                    print('SOCKET: Closed connection from: {}'.format(clients[notified_socket]['data'].decode('utf-8')))

                    # Remove from list for socket.socket()
                    sockets_list.remove(notified_socket)

                    # Remove from our list of users
                    del clients[notified_socket]

                    continue

                user = clients[notified_socket]
                print(f'SOCKET: Received message from {user["data"].decode("utf-8")}: {message["data"].decode("utf-8")}')

                # Send messages to clients
                for client_socket in clients:
                    if client_socket != notified_socket:
                        client_socket.send(id['header'] + id['data'] + message['header'] + message['data'])

