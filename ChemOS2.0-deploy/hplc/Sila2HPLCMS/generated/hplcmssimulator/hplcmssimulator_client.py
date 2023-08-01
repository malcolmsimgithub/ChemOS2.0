# Generated by sila2.code_generator; sila2.__version__: 0.10.3
# -----
# This class does not do anything useful at runtime. Its only purpose is to provide type annotations.
# Since sphinx does not support .pyi files (yet?), this is a .py file.
# -----

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Iterable, Optional

    from hplcmssimulator_types import (
        GetResults1st_IntermediateResponses,
        GetResults1st_Responses,
        Status_Responses,
        SubmitJob_IntermediateResponses,
        SubmitJob_Responses,
        ValveStatus_Responses,
    )
    from sila2.client import ClientMetadataInstance, ClientObservableCommandInstanceWithIntermediateResponses


class HPLCMSsimulatorClient:
    """

    Runs the ChemSpeed platform using a given synthesis procedure

    """

    def Status(self, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None) -> Status_Responses:
        """
        Gets the Status of the HPLC
        """
        ...

    def ValveStatus(
        self, Purpose: str, StatusUpdate: str, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None
    ) -> ValveStatus_Responses:
        """
        Gets the Status of the HPLC
        """
        ...

    def SubmitJob(
        self, JobFile: str, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None
    ) -> ClientObservableCommandInstanceWithIntermediateResponses[SubmitJob_IntermediateResponses, SubmitJob_Responses]:
        """
        Submits a Job for the HPLCMS to do
        """
        ...

    def GetResults1st(
        self, JobFile: str, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None
    ) -> ClientObservableCommandInstanceWithIntermediateResponses[
        GetResults1st_IntermediateResponses, GetResults1st_Responses
    ]:
        """
        Gets data from characterization_1st job
        """
        ...
