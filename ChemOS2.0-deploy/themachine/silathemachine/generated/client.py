# Generated by sila2.code_generator; sila2.__version__: 0.10.1
from __future__ import annotations

from typing import Set

from sila2.client import SilaClient
from sila2.framework import FullyQualifiedFeatureIdentifier

from .themachine import TheMachineClient


class Client(SilaClient):

    TheMachine: TheMachineClient

    _expected_features: Set[FullyQualifiedFeatureIdentifier] = {
        FullyQualifiedFeatureIdentifier("org.silastandard/core/SiLAService/v1"),
        FullyQualifiedFeatureIdentifier("org.silastandard/examples/TheMachine/v1"),
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
