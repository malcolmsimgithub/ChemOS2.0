from .ChemSpeedModules import ChemSpeedModule
from .HardwareExceptions import *
from pathlib import Path
from typing import Iterable


class FiltrationRack(ChemSpeedModule):
    """
    FiltrationRack Object
    Describes the positions on the SPE / filtration rack (status and availability).
    Possible states: draw, collect, waste
    """
    states: dict = {
        "draw": "_D",
        "collect": "_C",
        "waste": "_W"
    }

    def __init__(self, name: str, role: str, output_path: Path, positions: Iterable, output_folder: str, defaults_path: Path, **kwargs):
        """
        Instantiates the FiltrationRack object by calling the __init__ method of the ChemSpeedModule MetaClass.
        Overrides the self.positions parameter (integer keys).
        """
        super().__init__(name, role, output_path, positions, output_folder, defaults_path)
        self._load_status()

    def get_filtration_position(self, identifier: str) -> int:
        """
        Returns an available filtration / SPE position.
        Sets the value in self.position to the corresponding experiment identifier.
        Raises a PositionNotAvailableError if no position is availabe.
        """
        self._load_status()
        for position in self.positions:
            if self.positions[position] is None:
                self.positions[position] = identifier
                self._save_status()
                return position

        raise PositionNotAvailableError()

    def get_full_position(self, vial_number: int, state: str = "draw") -> str:
        """
        Get a full zone name for a position on the SPE Rack.
        """
        return f"{self.name}{self.states[state]}:{vial_number}"

    def update_from_file(self) -> None:
        """
        Calls the private method of updating from the specified json file.
        """
        self._load_status()
