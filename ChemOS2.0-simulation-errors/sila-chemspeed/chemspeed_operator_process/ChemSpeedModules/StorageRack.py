from .ChemSpeedModules import ChemSpeedModule
from .HardwareExceptions import *
from ..Utils import load_json
from pathlib import Path
from typing import Union, Iterable
import os


class StorageRack(ChemSpeedModule):
    """
    Class to describe a storage rack on the ChemSpeed deck.
    Accounts for both solid (GDU cartridges) and liquid (20 mL, 8 mL, 4 mL) storage racks.
    """
    def __init__(self, name: str, role: str, output_path: Path, positions: Iterable, output_folder: str, defaults_path: Path, **kwargs):
        """
        Instantiates a storage rack object by calling the __init__ method of the ChemSpeedModule MetaClass.
        Sets the following additional attributes:
            - self.phase ("solid" or "liquid")
            - self.max_content
            - self.variable_pos (list): List of all vessel names that do not have a default compound assigned.
        """
        super().__init__(name, role, output_path, positions, output_folder, defaults_path)
        self.phase: str = kwargs["phase"]
        self.max_content = kwargs["capacity"]

        self.variable_pos = []
        self.restore_default()

    def _load_permanent_compounds(self):
        """
        Loads the permanent compounds on the rack from the respective json file in the default settings folder.
        Sets the respective values in self.positions.
        """
        permanent_compounds = load_json(os.path.join(self.defaults_path,"Permanent_Compounds.json"))[self.phase]
        for compound, position in zip(permanent_compounds, permanent_compounds.values()):
            if position in self.positions:
                self.positions[position] = compound

    def restore_default(self) -> None:
        """
        Loads the default rack configuration from self.default_file.
        Resets the list of variable positions.
        """
        self.positions = {position: None for position in self.positions}
        self._load_permanent_compounds()
        self.variable_pos = [position for position in self.positions if self.positions[position] is None]
        self._save_status()

    def update_storage(self, position: str, compound: str) -> None:
        """
        Updates the self.positions dictionary based on incoming position-compound mapping.
        """
        if position not in self.positions:
            raise PositionNotAvailableError()

        self.positions[position] = compound
        self._save_status()

    def get_positions(self) -> list:
        """
        Returns a list of all positions on the rack.
        """
        return list(self.positions.keys())

    def get_position_content(self, vessel: str) -> Union[str, None]:
        """
        Returns the value for a position key on the storage rack.
        """
        if vessel not in self.positions:
            raise PositionNotAvailableError()
        return self.positions[vessel]

    def get_variable_positions(self) -> list:
        """
        Returns a list of all variable positions on the rack (i.e. the attribute self.variable_pos).
        """
        return self.variable_pos
