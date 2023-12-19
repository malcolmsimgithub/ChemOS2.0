# This file is adapted from a paramiko demo, and thus licensed under LGPL 2.1.
# Original Copyright (C) 2003-2007  Robey Pointer <robeypointer@gmail.com>
# Edits Copyright (C) 2010 The IPython Team
# Edits Copyright (C) 2020 Theophile Gaudin
#
# Paramiko is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Paramiko is distrubuted in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Paramiko; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02111-1301  USA.
from paramiko import SSHClient, ProxyCommand
import paramiko
import logging
import socketserver
import select
import socket


class ForwardServer (socketserver.ThreadingTCPServer):
    daemon_threads = True
    allow_reuse_address = True


class Handler(socketserver.BaseRequestHandler):
    logger = logging.getLogger(__name__)

    def handle(self):
        try:
            channel = self.ssh_transport.open_channel('direct-tcpip',
                                                   (self.chain_host, self.chain_port),
                                                   self.request.getpeername())
        except Exception as e:
            self.logger.debug('Incoming request to %s:%d failed: %s' % (self.chain_host,
                                                              self.chain_port,
                                                              repr(e)))
            return

        if channel is None:
            self.logger.debug('Incoming request to %s:%d was rejected by the SSH server.' %
                    (self.chain_host, self.chain_port))
            return

        self.logger.debug('Connected!  Tunnel open %r -> %r -> %r' % (self.request.getpeername(),
                                                            channel.getpeername(), (self.chain_host, self.chain_port)))
        while True:
            r, w, x = select.select([self.request, channel], [], [])
            if self.request in r:
                data = self.request.recv(1024)
                if len(data) == 0:
                    break
                channel.send(data)
            if channel in r:
                data = channel.recv(1024)
                if len(data) == 0:
                    break
                self.request.send(data)
        channel.close()
        self.request.close()
        self.logger.debug('Tunnel closed')


def tunnel_factory(local_port, remote_host, remote_port, transport):
    # this is a little convoluted, but lets me configure things for the Handler
    # object.  (SocketServer doesn't give Handlers any way to access the outer
    # server normally.)
    class SubHandler(Handler):
        chain_host = remote_host
        chain_port = remote_port
        ssh_transport = transport
    return ForwardServer(('localhost', local_port), SubHandler)


def connect(username, hostname, key_filename, proxycommand=None):
    ssh = SSHClient()
    ssh.set_missing_host_key_policy(paramiko.client.AutoAddPolicy)
    proxy = ProxyCommand(proxycommand) if proxycommand is not None else None
    ssh.connect(username=username,
                hostname=hostname,
                key_filename=key_filename,
                look_for_keys=False,
                sock=proxy)
    return ssh.get_transport()


def get_free_tcp_port():
    # TODO this is quite hacky, it has race condition problem
    tcp = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    tcp.bind(('', 0))
    addr, port = tcp.getsockname()
    tcp.close()
    return port
