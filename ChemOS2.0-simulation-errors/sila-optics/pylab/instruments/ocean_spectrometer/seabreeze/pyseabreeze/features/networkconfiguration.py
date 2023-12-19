from seabreeze.pyseabreeze.features._base import SeaBreezeFeature


# Definition
# ==========
#
# TODO: This feature needs to be implemented for pyseabreeze
#
class SeaBreezeNetworkConfigurationFeature(SeaBreezeFeature):
    identifier = "network_configuration"

    def get_number_of_network_interfaces(self):
        raise NotImplementedError("implement in derived class")

    def get_network_interface_connection_type(self, interface_index):
        raise NotImplementedError("implement in derived class")

    def get_network_interface_enable_state(self, interface_index):
        raise NotImplementedError("implement in derived class")

    def set_network_interface_enable_state(self, interface_index, enable_state):
        raise NotImplementedError("implement in derived class")

    def run_network_interface_self_test(self, interface_index):
        raise NotImplementedError("implement in derived class")

    def save_network_interface_connection_settings(self, interface_index):
        raise NotImplementedError("implement in derived class")
