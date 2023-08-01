from typing import Union, Any, Tuple
from ..Utils import timestamp_datetime

class Operation(object):
    """
    Class to describe an operation to be performed by the ChemSpeed.
    Contains:
        - the type of operation to be performed
        - operation-specific keyword arguments
    """
    available_operations: list = [
        "transfer_compound",
        "schlenk_cycle",
        "reflux",
        "prime_pumps",
        "request_confirmation",
        "unmount_all",
        "communicate",
        "set_drawer",
        "filter_collect",
        "submit_hplc_job",
        "inject_to_hplc"
    ]

    def __init__(
            self,
            experiment_identifier: list,
            task: str,
            parameters: dict = {},
            target: list = [],
            source: Union[str, None] = None,
            executed: Union[bool, str] = False
    ):
        """
        Instantiates an operation by setting the attributes self.task and self.parameters.
        Asserts that the task is within the scope of available tasks.
        """
        assert task in self.available_operations
        self.task: str = task
        self.parameters: dict = parameters
        self.experiment_identifier = experiment_identifier
        self.target = target
        self.source = source
        self.executed = executed

    def update_keyword_argument(self, keywords: dict) -> None:
        """
        Updates keyword arguments in self.parameters.
        Adds the keywords if not present.
        """
        for key in keywords:
            self.parameters[key] = keywords[key]

    def add_to_keyword_argument(self, keywords: dict) -> None:
        """
        Adds a value to a keyword argument.
        Only works if the keyword argument is either a list or a set.
        """
        for key in keywords:
            if type(self.parameters[key]) is list:
                self.parameters[key].append(keywords[key])
            if type(self.parameters[key]) is set:
                self.parameters[key].add(keywords[key])

    def group_operations(self, identifier: list, target: list) -> None:
        """
        Adds the parameters of a second, similar step (identifier & target_zone) to group them into one step
        """
        self.experiment_identifier = [*self.experiment_identifier, *identifier]
        self.target = [*self.target, *target]

    def is_similar(self, reference: Any, key: str) -> bool:
        """
        Compares the operation with another operation based on the task and a given keyword.
        Returns True if task and keyword argument are identical, otherwise returns False
        """
        if self.task != reference.task:
            return False
        if self.parameters[key] != reference.parameters[key]:
            return False

        return True

    def execute(self) -> Tuple[str, list, str, dict]:
        """
        Returns the task and the parameters to perform the execution of the operation.
        Sets the timestamp for execution in self.executed.
        """
        self._set_executed()
        return self.task, self.target, self.source, self.parameters

    def _set_executed(self) -> None:
        """
        Sets the self.executed variable to the current timestamp.
        """
        self.executed = timestamp_datetime()

    def get_dictionary(self) -> dict:
        """
        Returns a dictionary of all attributes.
        """
        operation_dict = {
            "task": self.task,
            "parameters": self.parameters,
            "experiment_identifier": self.experiment_identifier,
            "executed": self.executed
        }
        return operation_dict

    def is_in(self, attribute: str, comparison: Union[list, set]) -> bool:
        """
        Checks if any of the values within the comparison list are in the given attribute.
            - "experiment_identifier": self.experiment_identifier
            - "target_zone": self.target_zone
        """
        attributes = {
            "experiment_identifier": self.experiment_identifier,
            "target_zone": self.target
        }
        assert attribute in attributes
        for value in comparison:
            if value in attributes[attribute]:
                return True
        return False

    def get_compound(self) -> Union[dict, None]:
        """
        Returns the compound to be transferred in the current operation, if present.
        Returns None otherwise.
        """
        if "compound" in self.parameters:
            return {key: self.parameters[key] for key in ("compound", "specification")}
        else:
            return None


class Communication(Operation):
    """
    Specific operation to instantiate a communication operation object.
    """
    def __init__(self, message: str, thread: str):
        super().__init__(
            experiment_identifier=[None],
            task="communicate",
            parameters={"message": message, "thread": thread}
        )
