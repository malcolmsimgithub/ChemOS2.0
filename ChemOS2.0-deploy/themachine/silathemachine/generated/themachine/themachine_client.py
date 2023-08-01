# Generated by sila2.code_generator; sila2.__version__: 0.10.1
# -----
# This class does not do anything useful at runtime. Its only purpose is to provide type annotations.
# Since sphinx does not support .pyi files (yet?), so this is a .py file.
# -----

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:

    from typing import Iterable, Optional

    from sila2.client import ClientMetadataInstance, ClientObservableCommandInstanceWithIntermediateResponses
    from themachine_types import Runjob_IntermediateResponses, Runjob_Responses


class TheMachineClient:
    """

    TheMachine RemoteControl via Sila

    """

    def Runjob(
        self, Jobfile: str, ProcedureScript: str, *, metadata: Optional[Iterable[ClientMetadataInstance]] = None
    ) -> ClientObservableCommandInstanceWithIntermediateResponses[Runjob_IntermediateResponses, Runjob_Responses]:
        """
        Runs a job on the machine platform
        """
        ...
