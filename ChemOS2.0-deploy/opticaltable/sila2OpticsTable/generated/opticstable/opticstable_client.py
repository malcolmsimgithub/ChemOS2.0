# Generated by sila2.code_generator; sila2.__version__: 0.10.3
# -----
# This class does not do anything useful at runtime. Its only purpose is to provide type annotations.
# Since sphinx does not support .pyi files (yet?), this is a .py file.
# -----

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Iterable, Optional

    from opticstable_types import Status_Responses, SubmitJob_IntermediateResponses, SubmitJob_Responses
    from sila2.client import ClientMetadataInstance, ClientObservableCommandInstanceWithIntermediateResponses


class OpticsTableClient:
    """
    Simulates the optics table platform
    """

    def Status(self, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None) -> Status_Responses:
        """
        Gets the Status of the optics table
        """
        ...

    def SubmitJob(
        self, JobFile: str, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None
    ) -> ClientObservableCommandInstanceWithIntermediateResponses[SubmitJob_IntermediateResponses, SubmitJob_Responses]:
        """
        Submits a Job for the OpticsTableSimulator to do
        """
        ...
