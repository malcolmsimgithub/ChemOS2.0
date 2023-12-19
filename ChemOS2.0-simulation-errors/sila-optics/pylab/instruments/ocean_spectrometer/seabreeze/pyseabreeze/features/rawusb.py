from seabreeze.pyseabreeze.features._base import SeaBreezeFeature
from seabreeze.pyseabreeze.protocol import ProtocolInterface
from seabreeze.pyseabreeze.transport import USBTransport


class SeaBreezeRawUSBBusAccessFeature(SeaBreezeFeature):
    r"""
    Reading and writing raw usb

    Example usage

    .. code-block:: python

     >>> import struct  # needed for packing binary data into bytestrings
     >>> from seabreeze.spectrometers import Spectrometer
     >>> spec = Spectrometer.from_first_available()
     # features need to be accessed on the spectrometer instance via spec.f or spec.features
     # the data you provide needs to packed into a bytestring via struct.pack
     >>> spec.f.raw_usb_bus_access.raw_usb_write(data=struct.pack('<B', 0xFE), endpoint='primary_out')

     # when reading via the raw_usb feature you can easily block the spectrometer.
     # so make sure to only read when you expect new data (after you sent a command)
     # and only read as many bytes as the command returns
     >>> output = spec.f.raw_usb_bus_access.raw_usb_read(endpoint='primary_in', buffer_length=16)

     >>> print(output)
     b'\x00\x08p\x17\x00\x00\x00\x00\x00\x08\x00\x00\x00\x00\x80U'
     # extract by providing the data format from the developer datasheet and `output` to struct.unpack

    """
    identifier = "raw_usb_bus_access"
    _required_protocol_cls = ProtocolInterface

    @classmethod
    def supports_protocol(cls, protocol):
        return isinstance(protocol.transport, USBTransport)

    def raw_usb_write(self, data, endpoint):
        if endpoint not in {"primary_out", "secondary_out"}:
            raise ValueError("endpoint has to be one of 'primary_out', 'secondary_out'")
        if endpoint == "secondary_out":
            pass
        self.protocol.transport.write(data)

    def raw_usb_read(self, endpoint, buffer_length=None):
        if endpoint == "primary_in":
            return self.protocol.transport.read(size=buffer_length, mode="low_speed")
        elif endpoint == "secondary_in":
            return self.protocol.transport.read(size=buffer_length, mode="high_speed")
        elif endpoint == "secondary_in2":
            return self.protocol.transport.read(
                size=buffer_length, mode="high_speed_alt"
            )
        else:
            raise ValueError(
                "endpoint has to be one of 'primary_in', 'secondary_in', 'secondary_in2'"
            )
