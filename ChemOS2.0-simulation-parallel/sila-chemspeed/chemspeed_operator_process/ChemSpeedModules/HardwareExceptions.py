class StorageException(Exception):
    """
    Parent class for all exceptions related to compound storage.
    """
    def __init__(self, *args):
        super().__init__(*args)


class HardwareException(Exception):
    """
    Parent Class for all exceptions related to ChemSpeed hardware.
    """
    def __init__(self, *args):
        super().__init__(*args)


class UnknownHardwareException(Exception):
    """
    Error to be raised when an unknown hardware component is loaded.
    """
    def __init__(self, *args):
        super().__init__(*args)


class StorageFullError(StorageException):
    """
    Error to be raised when a specific storage rack has no free positions available.
    """
    def __init__(self, *args):
        super().__init__(*args)


class CompoundNotAvailableError(StorageException):
    """
    Error to be raised when a specific compound is not available from a storage.
    """
    def __init__(self, *args):
        super().__init__(*args)


class PositionNotAvailableError(StorageException):
    """
    Error to be raised when a specific storage position is not available.
    """
    def __init__(self, *args):
        super().__init__(*args)
