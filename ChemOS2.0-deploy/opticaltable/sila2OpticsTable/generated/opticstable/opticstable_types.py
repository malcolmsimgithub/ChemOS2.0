# Generated by sila2.code_generator; sila2.__version__: 0.10.3
from __future__ import annotations

from typing import NamedTuple


class Status_Responses(NamedTuple):
    Termination: str
    """
    Termination message
    """


class SubmitJob_Responses(NamedTuple):
    Termination: str
    """
    Termination message
    """


class SubmitJob_IntermediateResponses(NamedTuple):
    Status: str
    """
    Status of batch
    """

    Payload: bytes
    """
    Data Payload intermediate response
    """
