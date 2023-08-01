from pathlib import Path
from typing import List, Tuple

from ..Utils import load_json
from .Operation import Operation, Communication
from ..ChemSpeedModules import ISynth

import os


class StepManager(object):

    def __init__(self, default_settings: Path):
        self.step_settings: dict = load_json(os.path.join(default_settings,"Operations_Configuration.json"))

    def global_pre_operation(self, task: str, specification: str, thread: str) -> List[Operation]:
        """
        Returns a list of threads for the global pre-operation steps to be performed before the actual operation.
        """
        communication = Communication(f"Now performing the task {task}", thread)
        pre_op = [Operation([None], i["name"], i["parameters"]) for i in self.step_settings[task][specification]["operation"]]
        return [communication, *pre_op]

    def process_operations(self, operations: List[Operation], reactor: ISynth, group_by_reactor: bool = True) -> List[Operation]:
        """
        Processes a group of threads:
            - groups them by the reactor configuration, if desired
            - gets the pre-and post-operation steps to modify the drawer settings, if necessary
            - groups transfer steps within the threads
        """
        all_operations = []
        if group_by_reactor:
            groups = self._group_by_drawer(operations, reactor)
        else:
            groups = [operations]

        for group in groups:
            pre_operation, post_operation = self._get_drawer_operations(group, reactor)
            operations_processed = self._group_addition_steps(group)
            all_operations = all_operations + [*pre_operation, *operations_processed, *post_operation]

        return all_operations

    @staticmethod
    def _group_addition_steps(operations: List[Operation]) -> List[Operation]:
        """
        Groups addition steps (task = "transfer_compound") if the compounds transferred are identical.
        Avoids unnecessary wash steps and makes compound transfer more efficient.
        """
        operations_processed = []

        for operation in operations:
            try:
                for operation_proc in operations_processed:
                    if operation.is_similar(operation_proc, "compound"):
                        operation_proc.group_operations(operation.experiment_identifier, operation.target)
                        raise StopIteration
            except StopIteration:
                continue

            operations_processed.append(operation)

        return operations_processed

    def _get_drawer_operations(self, operations: List[Operation], reactor: ISynth) -> Tuple[list, list]:
        """
        Checks if setting the drawers is necessary to perform the actual operation (i.e. if settings are specified).
        Returns the operation to set the reactor drawer(s) into the correct position(s)
        before executing a list of threads.
        """
        target_settings = self.step_settings[operations[0].task][operations[0].parameters["specification"]]["drawers"]
        if not target_settings:
            return [], []

        identifiers = []
        for operation in operations:
            identifiers = [*identifiers, *operation.experiment_identifier]
        drawer_ids = reactor.get_active_drawers(identifiers)

        parameters = {"target_setting": target_settings}

        pre_operation = [Operation(identifiers, "set_drawer", parameters=parameters, target=list(drawer_ids))]
        post_operation = [Operation(identifiers, "set_drawer", parameters={}, target=list(drawer_ids))]

        return pre_operation, post_operation

    @staticmethod
    def _group_by_drawer(operations: List[Operation], reactor: ISynth) -> List[List[Operation]]:
        """
        ChemSpeed / ISynth-specific method
        Groups threads by the drawer which is part of the transfer step.

        ATTENTION: FAILS IF A STEP INVOLVES MORE THAN ONE VIAL ON THE ISYNTH (i.e. transfer within ISYNTH)!
        # TODO: Implement exception to account for that case.
        """
        grouped_operations = []
        drawers = reactor.get_all_drawers(reaction_ids=True)

        for drawer in drawers:
            drawer_operations = [operation for operation in operations if operation.is_in("experiment_identifier", drawer)]  # TODO: Change existence in drawer based on operation.target_zone
            operations = list(set(operations) - set(drawer_operations))
            grouped_operations.append(drawer_operations)

        return grouped_operations
