# Generated by sila2.code_generator; sila2.__version__: 0.10.1
from __future__ import annotations

from typing import TYPE_CHECKING

from sila2.server import MetadataDict, ObservableCommandInstanceWithIntermediateResponses

from ..generated.hplcmssimulator import (
    HPLCMSsimulatorBase,
    Status_Responses,
    SubmitJobAutosampler_IntermediateResponses,
    SubmitJobAutosampler_Responses,
    SubmitJobChemspeed_IntermediateResponses,
    SubmitJobChemspeed_Responses,
    ValveStatus_Responses,
)

if TYPE_CHECKING:
    from ..server import Server


class HPLCMSsimulatorImpl(HPLCMSsimulatorBase):
    def __init__(self, parent_server: Server) -> None:
        super().__init__(parent_server=parent_server)

    def Status(self, *, metadata: MetadataDict) -> Status_Responses:
        raise NotImplementedError  # TODO

    def ValveStatus(self, Purpose: str, StatusUpdate: str, *, metadata: MetadataDict) -> ValveStatus_Responses:
        raise NotImplementedError  # TODO

    def SubmitJobAutosampler(
        self,
        JobFile: str,
        *,
        metadata: MetadataDict,
        instance: ObservableCommandInstanceWithIntermediateResponses[SubmitJobAutosampler_IntermediateResponses],
    ) -> SubmitJobAutosampler_Responses:
        instance.begin_execution()  # set execution status from `waiting` to `running`
        raise NotImplementedError  # TODO

    def SubmitJobChemspeed(
        self,
        JobFile: str,
        *,
        metadata: MetadataDict,
        instance: ObservableCommandInstanceWithIntermediateResponses[SubmitJobChemspeed_IntermediateResponses],
    ) -> SubmitJobChemspeed_Responses:
        instance.begin_execution()  # set execution status from `waiting` to `running`
        raise NotImplementedError  # TODO