import copy
from ..Utils import load_json
from .Operation import Operation
from pathlib import Path
from typing import Any, Union, Dict, List


class Synthesis(object):
    """
    Class to describe a specific synthesis procedure.
    """

    def __init__(self, procedure_file: Path, thread="Synthesis", placeholder_pattern: str = "$"):
        self.steps = dict()
        self.compounds = {"solid": [], "liquid": []}
        self.variables = set()
        self.placeholder_pattern = placeholder_pattern
        self.thread = thread
        self._parse_procedure(procedure_file)

    def _parse_procedure(self, procedure_file: Path) -> None:
        """
        Loads the synthesis procedure from the specified json file.
        """
        procedure = load_json(procedure_file)

        for step in procedure:
            variables = self._parse_step(procedure[step])
            self.steps[step] = {
                "task": procedure[step]["task"],
                "local": procedure[step]["local"],
                "pre-operation": procedure[step]["pre-operation"],
                "variable": variables,
                "parameters": procedure[step]["parameters"]
            }

    def _check_for_placeholder(self, variable: Any):
        """
        Checks if a specific variable (dictionary key) is a placeholder:
            - type must be str
            - must contain placeholder pattern
        """
        if isinstance(variable, str) and self.placeholder_pattern in variable:
            return True
        else:
            return False

    def _parse_step(self, reaction_step: dict) -> dict:
        """
        Parses a single reaction step by iterating through the parameters and specifying variables.
        Returns the operation variables as a dictionary of dictionaries: e.g.
        e.g.variables = {
                $A$: {"type": "compound", "specification": "solid"},
                $vol$: {"type": "quantity", "specification": "liquid"}
            }.
        """
        variables = dict()

        for parameter in reaction_step["parameters"]:
            if self._check_for_placeholder(reaction_step["parameters"][parameter]):
                name = reaction_step["parameters"][parameter]
                variables[name] = {
                    "type": parameter,
                    "specification": reaction_step["parameters"]["specification"]
                }

        self._add_variables(variables)
        return variables

    def _add_variables(self, variables: dict) -> None:
        """
        Adds all variables to self.variables.
        Adds all variables of the type "compound" to the self.compounds dictionary.
        """
        for variable, variable_props in zip(variables, variables.values()):
            self.variables.add(variable)
            if variable_props["type"] == "compound":
                self.compounds[variable_props["specification"]].append(variable)
        # TODO: Double-check where the self.compounds dictionary is actually necessary

    def get_compound_types(self) -> list:
        """
        Returns a list of all placeholder compound types found in the current procedure.
        """
        compound_types = []
        for compounds in self.compounds.values():
            for compound in compounds:
                compound_types.append(compound)
        return compound_types

    def get_operations(self, step_id: str, reaction_details: Dict[str, dict]) -> List[Operation]:
        """
        Gets a list of specific threads for a given reaction step and a list of reactions.
        Replaces the placeholder variables by actual compound identifiers.
        Returns a list of operation objects.
        """
        if self.step_is_local(step_id):
            operations = [self._generate_operation(self.steps[step_id], key, reaction_details[key]) for key in reaction_details]
        else:
            operations = [self._generate_operation(self.steps[step_id], None, {"vial": None})]

        return operations

    @staticmethod
    def _generate_operation(step: dict, identifier: Union[str, None], mapping: dict) -> Operation:
        """
        Private method which converts an abstract step into an actual operation for a specific reaction.
        Substitutes the placeholders in step["parameters"] according to the mapping provided
        Returns the Operation object.
        """
        step = copy.deepcopy(step)
        for variable in step["variable"]:
            var_type = step["variable"][variable]["type"]
            step["parameters"][var_type] = mapping[variable]

        return Operation([identifier], step["task"], step["parameters"], [mapping["vial"]])

    def step_is_local(self, step_id: str) -> bool:
        """
        Checks if a step from the loaded procedure is
            - local (individual operation for each reaction, e.g. compound transfer) -> True
            - global (same operation for all reactions, e.g. heat / stir) -> False
        """
        return self.steps[step_id]["local"]

    def get_step_details(self, step_id: str) -> dict:
        """
        Gets the details (taks and specification) for a specific synthesis step
        to be able to clearly address it in the step settings.
        """
        details = {
            "task": self.steps[step_id]["task"],
            "specification": self.steps[step_id]["parameters"]["specification"],
            "thread": self.thread
        }
        return details

    def perform_global_pre_operation(self, step_id) -> bool:
        """
        Checks if a global pre-operation (e.g. pump priming) is to be performed for a given step.
        """
        return self.steps[step_id]["pre-operation"]
