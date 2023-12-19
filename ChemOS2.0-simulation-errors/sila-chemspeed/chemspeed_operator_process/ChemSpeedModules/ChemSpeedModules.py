from abc import ABCMeta
from pathlib import Path
from typing import Iterable
from ..Utils import load_json, save_json
import os


class ChemSpeedModule(metaclass=ABCMeta):
    """
    Abstract Base Class for a ChemSpeed Module object.
    """
    def __init__(self, name: str, role: str, output_path: Path, positions: Iterable, output_folder: str, defaults_path: Path):
        self.name = name
        self.role = role
        self.output_path = os.path.join(output_path,output_folder)
        self.output_file = os.path.join(self.output_path,f"{self.name}.json")
        self.defaults_path = defaults_path

        self.positions = {f"{self.name}:{i}": None for i in positions}

    def _load_status(self) -> None:
        """
        Updates self.positions by all data specified in the output file.
        If the output file does not exist, dictionary update is skipped.
        """
        try:
            loaded_positions = load_json(self.output_file)
            self.positions.update(loaded_positions)
        except FileNotFoundError:
            self._save_status()

    def _save_status(self) -> None:
        """
        Saves the current status of the module into a json file.
        """
        save_json(self.positions, self.output_file)
