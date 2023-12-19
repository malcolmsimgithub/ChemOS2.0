from pathlib import Path
from typing import Iterable

from .FiltrationRack import FiltrationRack
from .HardwareExceptions import *
from .HPLCInterface import HPLCInterface
from .ISynth import ISynth
from .StorageRack import StorageRack


class ChemSpeedFactory(object):
    """
    Factory Pattern Class for generating ChemSpeed Module Objects.
    """
    def __init__(self):
        self._builders = {
            "FiltrationRack": self._get_filtration_rack,
            "ISynth": self._get_isynth,
            "HPLCInterface": self._get_hplc_interface,
            "StorageRack": self._get_storage_rack,
        }

    def __call__(self, object_type: str, **kwargs):
        return self.get(object_type, **kwargs)

    def get(
            self,
            object_type: str,
            name: str,
            role: str,
            output_path: Path,
            positions: str,
            output_folder: str,
            defaults_path: Path,
            **kwargs
    ):
        """
        Product of the factory pattern to generate the ChemSpeed objects.
        """
        if object_type not in self._builders:
            raise UnknownHardwareException(object_type)

        return self._builders[object_type](
            name=name,
            role=role,
            output_path=output_path,
            positions=eval(positions),
            output_folder=output_folder,
            defaults_path=defaults_path,
            **kwargs
        )

    @staticmethod
    def _get_filtration_rack(**kwargs):
        """
        Creator for the FiltrationRack class.
        """
        return FiltrationRack(**kwargs)

    @staticmethod
    def _get_hplc_interface(**kwargs):
        """
        Creator for the HPLCInterface class.
        """
        return HPLCInterface(**kwargs)

    @staticmethod
    def _get_isynth(**kwargs):
        """
        Creator for the ISynth class.
        """
        return ISynth(**kwargs)

    @staticmethod
    def _get_storage_rack(**kwargs):
        """
        Creator for the StorageRack class.
        """
        return StorageRack(**kwargs)
