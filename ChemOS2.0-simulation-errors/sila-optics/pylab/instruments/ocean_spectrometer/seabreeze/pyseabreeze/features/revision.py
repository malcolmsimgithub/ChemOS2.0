from seabreeze.pyseabreeze.features._base import SeaBreezeFeature


# Definition
# ==========
#
# TODO: This feature needs to be implemented for pyseabreeze
#
class SeaBreezeRevisionFeature(SeaBreezeFeature):
    identifier = "revision"

    def hardware_revision(self):
        raise NotImplementedError("implement in derived class")

    def revision_firmware_get(self):
        raise NotImplementedError("implement in derived class")
