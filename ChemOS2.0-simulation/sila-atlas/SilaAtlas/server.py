# Generated by sila2.code_generator; sila2.__version__: 0.10.1

from typing import Optional
from uuid import UUID

from sila2.server import SilaServer

from .feature_implementations.atlas_impl import AtlasImpl
from .generated.atlas import AtlasFeature


class Server(SilaServer):
    def __init__(self, server_uuid: Optional[UUID] = None):
        # TODO: fill in your server information
        super().__init__(
            server_name="TODO",
            server_type="TODO",
            server_version="0.1",
            server_description="TODO",
            server_vendor_url="https://gitlab.com/SiLA2/sila_python",
            server_uuid=server_uuid,
        )

        self.atlas = AtlasImpl(self)
        self.set_feature_implementation(AtlasFeature, self.atlas)
