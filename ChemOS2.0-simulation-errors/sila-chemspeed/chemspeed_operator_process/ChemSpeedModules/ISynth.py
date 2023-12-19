from pathlib import Path
import random
from typing import Union, Iterable, Set, List
from .ChemSpeedModules import ChemSpeedModule


class ISynth(ChemSpeedModule):
    """
    Class to describe the ISynth Reactor including all its reaction positions.
    Groups the positions by Drawers.
    """
    def __init__(self, name: str, role: str, output_path: Path, positions: Iterable, output_folder: str, defaults_path: Path, **kwargs):
        """
        Initializes the ISynth object by reading drawers and positions from ISynth.json.
        Sets the following attributes
            - self.drawers: list[dict]
            - self.positions: dict
        Information stored in both attributes is redundant, but useful to have.
        """
        super().__init__(name, role, output_path, positions, output_folder, defaults_path)

        self.drawers: list[dict] = [{} for _ in range(int(len(positions) / kwargs["drawer_size"]))]
        self._set_drawers(kwargs["drawer_size"])

    def set_batch_vials(self, num_rxn: int) -> None:
        """
        Specifies the positions to be used in a certain batch
        by setting the values to True in self.drawers and self.positions.

        Uses the first n positions by default to minimize the amount of active drawers.
        """
        for i, vial in enumerate(self.positions):
            if i < num_rxn:
                self.positions[vial] = True

        self._update_drawers()

    def _set_drawers(self, drawer_size: int) -> None:
        """
        Method for initialization of the ISynth object.
        Sets up self.drawers as a list of dictionaries of individual positions.
        !!! ATTENTION: REQUIRES DRAWER-WISE ORDER OF PASSED POSITIONS!!!
        """
        for i, position in enumerate(self.positions.keys()):
            drawer_index = int(i / drawer_size)
            self.drawers[drawer_index][position] = None

    def _update_drawers(self) -> None:
        """
        Updates the vial status in self.drawers based on the vial status in self.positions.
        """
        for drawer in self.drawers:
            for vial in drawer:
                drawer[vial] = self.positions[vial]

    def _update_vials(self) -> None:
        """
        Updates the vial status in self.positions based on the vial status in self.drawers.
        """
        self.positions = {}
        for drawer in self.drawers:
            self.positions = {**self.positions, **drawer}

    def _available_positions(self) -> list:
        """
        Returns the number of available (i.e. value is True) positions per drawer.
        """
        positions = [sum(value is True for value in drawer.values()) for drawer in self.drawers]
        return positions

    def assign_group(self, reaction_ids: set, assignments: Union[dict, None] = {}) -> dict:
        """
        Assigns a group of reaction identifiers to the corresponding positions using a minimum amount of drawers.
            - Checks if there is a drawer with enough available positions.
                If True, assigns the reactions.
                If False, splits the group of reactions according to availability and calls the function recursively.
            - Sets self.drawers and self.positions accordingly
        """
        positions = self._available_positions()
        for index, no_pos in enumerate(positions):
            if no_pos >= len(reaction_ids):
                assignments.update(self._assign_to_drawer(index, reaction_ids))
                return assignments

        max_available = max(positions)
        subset = set(random.sample(reaction_ids, max_available))
        assignments.update(self._assign_to_drawer(positions.index(max_available), subset))
        reaction_ids = reaction_ids - subset
        if reaction_ids:
            return self.assign_group(reaction_ids)

    def _assign_to_drawer(self, drawer_index: int, reaction_ids: set) -> dict:
        """
        Assigns a batch of reactions to the positions in the respective drawer.
        Assumes Pre-Assertion that the number of available positions is sufficient.
        Sets the values in self.drawers and self.positions
        """
        assignments = {}
        for reaction in reaction_ids:
            for vial in self.drawers[drawer_index]:
                if self.drawers[drawer_index][vial] is True:
                    self.drawers[drawer_index][vial] = reaction
                    assignments[reaction] = vial
                    break

        self._update_vials()
        return assignments

    def _get_vial(self, reaction_id) -> Union[str, None]:
        """
        Gets the vial identifier for a given reaction ID.
        Returns None if the reaction is not assigned to any vial.
        """
        for vial in self.positions:
            if self.positions[vial] == reaction_id:
                return vial

        return None

    def get_all_drawers(self, reaction_ids: bool) -> List[set]:
        """
        Returns a list of sets of vial identifiers or reaction identifiers (per drawer).
        """
        if reaction_ids:
            return [set(drawer.values()) for drawer in self.drawers if any(value is not None for value in drawer.values())]
        else:
            return [set(drawer.keys()) for drawer in self.drawers]

    def get_active_drawers(self, reaction_ids: Union[list, set]) -> Set[str]:
        """
        Returns a list of vial identifiers representing the active drawers containing the reactions given.
        """
        active_drawers = set()
        for reaction in reaction_ids:
            for drawer in self.drawers:
                if reaction in drawer.values():
                    active_drawers.add(list(drawer.keys())[0])

        return active_drawers

    def all_active_drawers(self) -> List[str]:
        """
        Returns a list of vial identifiers representing the active (i.e. reaction-containing) drawers.
        """
        active_drawers = []
        for drawer in self.drawers:
            if any(value is not None for value in drawer.values()):
                active_drawers.append(list(drawer.keys())[0])

        return active_drawers

    def clear_reactor(self):
        """
        Clears all assigned reactions from the reactor
        """
        for drawer in self.drawers:
            for position in drawer:
                drawer[position] = None
        self._update_vials()
