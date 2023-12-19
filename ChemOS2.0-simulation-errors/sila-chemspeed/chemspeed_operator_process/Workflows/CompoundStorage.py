from ..Utils import load_json, save_json
from pathlib import Path
from ..ChemSpeedModules.StorageRack import StorageRack
from ..ChemSpeedModules.HardwareExceptions import *
from ..Workflows import Operation
from typing import Callable, Union, List
import os
import pathlib


class CompoundStorage(object):
    """
    Class to handle compound storage and inventories
    """
    def __init__(self, components: dict, defaults_path: Path, output_path: Path, ui_confirmation: Callable):
        """
        Instantiates the CompoundStorage.

        Sets the following attributes:
            - self.storage_racks (list): list of dicts of StorageRack objects
            - self.ui_confirmation (Callable): Method to request user confirmations via the user interface.
            - self.output_file (Path): Path to the output file where the current inventory is saved.
            - self.storage (dict): Dictionary of compound names and current storage location
        """
        self.storage_racks = {"solid": {}, "liquid": {}}
        self._setup_storage_racks(components)

        self.current_inventory = {"solid": None, "liquid": None}
        self.inventories = {"solid": None, "liquid": None}

        self.ui_confirmation = ui_confirmation
        self.output_file = os.path.join(output_path,"Inventory.json")

        self.defaults_path = defaults_path

        self.storage = self._load_storage()
        self._update_storage_from_racks()

    def _setup_storage_racks(self, hardware: dict):
        """
        Iterates over all hardware modules, identifies StorageRacks and assigns them to self.storage_racks
        """
        for module in hardware.values():
            if type(module) is StorageRack and module.role == "storage":
                self.storage_racks[module.phase][module.name] = module

    def _save_storage(self) -> None:
        """
        Saves the self.storage dictionary to self.output_file.
        """
        save_json(self.storage, self.output_file)

    def _load_storage(self) -> dict:
        """
        Private method.
        Loads and returns the previously saved compound dictionary from the self.output_file.
        Returns an empty dictionary if the file does not exist.
        """
        if os.path.isfile(self.output_file):
            return load_json(self.output_file)
        else:
            return dict()

    def update_storage(self, **kwargs) -> None:
        """
        Updates self.storage based on all compounds passed as kwargs (name = {key: value, ...}
        """
        for compound in kwargs:
            self.update_compound(compound, **kwargs[compound])

    def _update_storage_from_racks(self):
        """
        Goes through all compounds in the racks and assigns the corresponding position in self.storage.
        If a compound is added to self.storage for the first time, the quantity is set to the rack max. quantity.
        """
        for phase in self.storage_racks:
            for rack in self.storage_racks[phase].values():
                for vessel in rack.get_positions():
                    compound = rack.get_position_content(vessel)
                    if compound in self.storage:
                        self.update_compound(compound, location=vessel)
                    elif compound:
                        self._create_compound(compound, location=vessel, phase=phase, quantity=rack.max_content)

    def get_compound_location(self, name: str) -> str:
        """
        Gets the name of the zone where a compound is located.

        If the compound has a specified location, the location is returned.
        Otherwise, the next inventory is set.
        """
        location = self.storage[name]["location"]
        phase = self.storage[name]["phase"]

        if location is not None:
            return location

        else:
            self.set_next_inventory(phase)
            return self.get_compound_location(name)

    def update_compound(self, name: str, **kwargs) -> None:
        """
        Updates the entry for the compound in self.storage from the passed keywords.
        Creates a new entry if it doesn't exist already.
        """
        if name not in self.storage:
            self._create_compound(name, **kwargs)
            return

        self.storage[name].update(kwargs)
        self._save_storage()

    def _create_compound(self, name: str, **kwargs) -> None:
        """
        Creates a new entry for a compound.
        """
        cmpd_entry = {
            "location": None,
            "phase": None,
            "quantity": 1000000
        }
        cmpd_entry.update(kwargs)
        self.storage[name] = cmpd_entry
        self._save_storage()

    def generate_inventories(self, operations: List[Operation]) -> None:
        """
        Generates a list of Inventory objects for each phase (solid, liquid) from a list of Operations
        that need to be executed sequentially.

        Each inventory object contains a static assignment of compounds and storage locations on the instrument deck.

        1. Extracts a list of all free, variable positions on all StorageRacks
        2. Sets up the list of inventories for each phase, as well as the current_inventory object.
        3. Iterates over all operations, extracts the compound (if any) and assigns it to the current inventory
           Once full, it is appended to all_inventories, and the next blank inventory is generated.
        4. Sets the iterators of all_inventories to self.inventories
        """
        permanent_compounds: dict = load_json(os.path.join(self.defaults_path,"Permanent_Compounds.json"))

        # get the variable (non-permanent) storage positions for all available racks for each phase
        variable_positions: dict = dict()
        for phase in self.storage_racks:
            var_pos_per_phase: list = []
            for rack in self.storage_racks[phase].values():
                var_pos_per_phase = [*var_pos_per_phase, *rack.get_variable_positions()]
            variable_positions[phase] = var_pos_per_phase

        # Set up the objects for all inventories to be generated, and the current inventory
        all_inventories = {phase: [] for phase in self.storage_racks}
        current_inventory = {phase: Inventory(variable_positions[phase]) for phase in self.storage_racks}

        # Extract compounds from each operation and add them to the current inventory until it is "full"
        for operation in operations:
            compound: Union[dict, None] = operation.get_compound()
            if compound:
                phase = compound["specification"]
                if not current_inventory[phase].contains(compound["compound"]) and compound["compound"] not in permanent_compounds[phase]:
                    self.update_compound(compound["compound"], phase=phase)
                    try:
                        current_inventory[phase].add_compound(compound["compound"])
                    except PositionNotAvailableError:
                        all_inventories[phase].append(current_inventory[phase])
                        current_inventory[phase] = Inventory(variable_positions[phase])
                        current_inventory[phase].add_compound(compound["compound"])

        for phase in current_inventory:
            all_inventories[phase].append(current_inventory[phase])
            self.inventories[phase] = iter(all_inventories[phase])
            self.set_next_inventory(phase)

    def set_next_inventory(self, phase):
        """
        Takes the next inventory from the self.inventories iterator for the corresponding phase.
        Sets the corresponding value in self.current_inventory.
        Modifies self.storage and self.storage_racks accordingly.
        """
        try:
            self.current_inventory[phase] = next(self.inventories[phase])
            for position in self.current_inventory[phase].all_positions():
                compound = self.current_inventory[phase].get_compound(position)
                self.update_compound(compound, location=position)

                for rack_name in self.storage_racks[phase]:
                    if rack_name in position:
                        self.storage_racks[phase][rack_name].update_storage(position, compound)

            self.ui_confirmation(message=f"The {phase} compound storage has been reset."
                                         f"Please Update the settings accordingly.")

        except StopIteration:
            self.current_inventory[phase] = None
            self.inventories[phase] = None
            raise ValueError("No matching inventory found. ")

    def clear_storage(self) -> None:
        """
        Resets all storage racks to default.
        Clears the current storage (no compound locations except for defaults).
        """
        for compound in self.storage:
            self.update_compound(compound, location=None)

        for phase in self.storage_racks:
            for rack in self.storage_racks[phase].values():
                rack.restore_default()

        self._update_storage_from_racks()


class Inventory(object):
    """
    Temporary assignment of compounds to reactors.
    Contains a dictionary mapping every reactor id to a compound name.
    """
    def __init__(self, positions: list):
        """
        Sets the attribute self.inventory as a dictionary mapping each position to a compound (default: None).
        """
        self.inventory = {position: None for position in positions}

    def add_compound(self, cmpd: str) -> None:
        """
        Adds a compound to the next free position in self.inventory.
        """
        free_position = self._get_next_free_position()
        self.inventory[free_position] = cmpd

    def all_positions(self) -> list:
        """
        Returns a list of all positions in the current inventory.
        """
        return list(self.inventory.keys())

    def get_compound(self, position: str) -> str:
        """
        Gets the compound assigned to a specific position in the inventory.
        """
        return self.inventory[position]

    def _get_next_free_position(self) -> str:
        """
        Iterates over self.inventory and returns the first position where the value is None.
        Raises a ValueError if no free position is available.
        """
        for position in self.inventory:
            if self.inventory[position] is None:
                return position

        raise PositionNotAvailableError

    def contains(self, compound: str) -> bool:
        """
        Checks whether a specific compound is in the values of self.inventory.
        """
        return compound in list(self.inventory.values())
