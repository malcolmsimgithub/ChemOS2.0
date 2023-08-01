
from .ChemSpeedModules import HPLCInterface, ISynth
from .Chemspeedcontroller.controller import ChemspeedController
import time
from pathlib import Path
from typing import Union, Callable
import os


class ChemSpeedExecutor(object):
    """
    Controller object for executing physical operations on the ChemSpeed Platform.
    Interface / higher-level API for the ChemspeedController.
    """
    def __init__(
            self,
            command_path: Path,
            output_path: Path,
            compound_location_mapping: Callable,
            hplc_controller: HPLCInterface,
            reactor: ISynth,
            simulation: bool = False
    ):
        """
        Parameters:
            command_path: Path to the ChemSpeed Folder
            output_path: Path where the log file is to be located.
            compound_location_mapping: Function to get
            simulation: Run Python-internal simulation
        """
        self.controller = ChemspeedController(command_path, logfile=os.path.join(output_path,"operation_log.log"), simulation=simulation)
        self.get_compound_location = compound_location_mapping
        self.hplc_interface = hplc_controller
        self.reactor = reactor
        self.status = "idle"
        self.simulation = simulation

    def _operate_if_idle(
            target_status: Union[str, None] = None
    ) -> Callable:
        """
        Decorator method to be used for any physical operation:
            - waits if instrument status is idle
            - sets instrument status to "busy"
            - if stirring: temporarily switches off stirring
            - performs physical operation(s) -> executes decorated function
            - if stirring: re-starts stirring
            - sets instrument status back to before (or to the target_status, if defined)

        IMPORTANT:  Only to be applied to functions that perform actual controller operations
                    (not to functions that call other functions that perform controller operations...)

        Parameters:
            target_status: Status to set the instrument after completing the operation.

        Returns:
            decorator
        """
        def decorator(function):
            def wrapper(self, *args, **kwargs):
                while self.status not in ("idle", "stirring"):
                    time.sleep(1)

                if self.status == "idle":
                    self.status = "busy"
                    function(self, *args, **kwargs)
                    while self.status == "interrupt":
                        time.sleep(1)
                    if target_status:
                        self.status = target_status
                    else:
                        self.status = "idle"

                elif self.status == "stirring":
                    self.status = "busy"
                    self.controller.set_stir(stir_zone="ISYNTH", state="off")
                    function(self, *args, **kwargs)
                    while self.status == "interrupt":
                        time.sleep(1)
                    if target_status == "idle":
                        self.status = target_status
                        self.controller.set_stir(stir_zone="ISYNTH", state="off")
                    elif target_status == "stirring":
                        self.status = target_status
                        self.controller.set_stir(stir_zone="ISYNTH", state="on", rpm=300)
                    else:
                        self.status = "stirring"
                        self.controller.set_stir(stir_zone="ISYNTH", state="on", rpm=300)
            return wrapper
        return decorator

    ####################################################
    # METHODS THAT DIRECTLY CORRESPOND TO AN OPERATION #
    ####################################################

    @_operate_if_idle(target_status=None)
    def transfer_compound(
            self,
            specification: str,
            target_zone: Union[list, str],
            quantity: float,
            source_zone: Union[str, None] = None,
            compound: Union[str, None] = None,
            dry: bool = False,
            **kwargs
    ) -> None:
        """
        Executes compound transfer on the instrument. Parses the required arguments for solid and liquid transfer.

        Args:
            specification: "solid" or "liquid"
            target_zone: Target zone on the instrument.
            quantity: Quantity of material transferred (in mL for liquids / mg for solids).
            compound (optional): Identifier of the compound to be transferred.
            source_zone (optional): Source zone of the compound to be transferred.
            dry (optional): True if the transfer should occur under water-free conditions
        """
        if not source_zone:
            source_zone = self.get_compound_location(compound)

        if specification == "solid":
            if kwargs["tolerance"]:
                chunk_size = float(kwargs["tolerance"])*2
            else:
                chunk_size = 1
            self.controller.transfer_solid(source=source_zone, destination=target_zone, weight=quantity, chunk=chunk_size)
        elif specification == "liquid":
            if dry:
                self.controller.transfer_liquid(
                    source=source_zone,
                    destination=target_zone,
                    volume=quantity,
                    src_flow=20,
                    dst_flow=40,
                    needle=1,
                    airgap=0.1,
                    post_airgap=0.1,
                    **kwargs
                )
            else:
                self.controller.transfer_liquid(
                    source=source_zone,
                    destination=target_zone,
                    volume=quantity,
                    src_flow=20,
                    dst_flow=40,
                    airgap=0.1,
                    post_airgap=0.1,
                    **kwargs
                )

    def schlenk_cycling(
            self,
            no_cycles: int,
            **kwargs
    ) -> None:
        """
        Performs Schlenk (evacuation to 1 mbar / backfill with inert gas) cycles on the ISynth.

        Args:
            no_cycles: Number of Schlenk cycles to perform.
        """
        drawers = self.reactor.all_active_drawers()
        for _ in range(no_cycles):
            self.set_drawer(target_zone=drawers, target_setting={"state": "close", "environment": "vacuum"})
            time.sleep(0.5)  # required, because otherwise communication between controller and manager can be too slow
            self.controller.set_isynth_vacuum(state="on", vacuum=1)
            self.controller.wait(60)
            self.set_drawer(target_zone=drawers, target_setting={"state": "close", "environment": "inert"})
            time.sleep(0.5)
            self.controller.set_isynth_vacuum(state="on", vacuum=1000)
            self.controller.wait(30)
            self.controller.set_isynth_vacuum(state="off")
        self.set_drawer(target_zone=drawers, target_setting={"state": "close", "environment": "none"})

    def reflux(
            self,
            temperature: float,
            heating: float,
            cooling: float,
            **kwargs
    ) -> None:
        """
        Sets ISynth under reflux (defined heating temperature, 0 °C as cooling temperature) for defined heating time,
        followed by defined cooling time.

        Args:
            temperature: Heating temperature (in °C).
            heating: Heating time (in h)
            cooling: Cooling time (in h).
        """
        self.set_reflux(temperature)
        self.controller.wait(duration=int(heating*3600))
        self.set_reflux(20)
        self.controller.wait(duration=int(cooling*3600))
        self.stop_reflux()

    @_operate_if_idle(target_status=None)
    def prime_pumps(
            self,
            **kwargs
    ) -> None:
        """
        Primes all pumps in order to avoid residual gas in the tubing.
        """
        for _ in range(2):
            self.controller.transfer_liquid(source="VALVEB:1;VALVEB:2;VALVEB:3;VALVEB:4", destination="WASTE1", volume=7.5, src_flow=20, dst_flow=40, rinse_volume=0)

    @_operate_if_idle(target_status=None)
    def set_drawer(
            self,
            target_zone: Union[str, list],
            target_setting: Union[dict, None] = None,
            **kwargs
    ) -> None:
        """
        Sets a specified ISYNTH drawer to the target_zone status.

        Args:
            target_zone: Zone describing the drawers to be set
            target_setting: Description of the target_zone setting, e.g. {"state": "closed", "environment": None}
        """
        if not target_setting:
            target_setting = {"state":  "close", "environment":  "none"}
        self.controller.set_drawer(zone=target_zone, **target_setting)

    def filter_collect(
            self,
            source_zone: str,
            target_zone: str,
            volume: float,
            pre_wash: Union[str, None] = None,
            pre_wash_volume: float = 0,
            post_wash: Union[str, None] = None,
            post_wash_volume: float = 0
    ) -> None:
        """
        Perform filtration of the solution on SPE Rack.

        Args:
            source_zone: Source zone of the solution to be filtered
            target_zone: Target zone on the SPE rack given as SPE:$NUMBER
            volume: Volume of liquid to be collected for filtration
            pre_wash: Name of the solvent for pre-washing
            pre_wash_volume: Volume of the solvent for pre-washing
            post_wash: Name of the solvent for post-washing
            post_wash_volume: Volume of the solvent for post-washing
        """
        # Quick, ugly fix for SPE:$NUMBER -> SPE_W:$NUMBER
        target_position: list = target_zone.split(":")
        waste_position = f"{target_position[0]}_W:{target_position[1]}"
        collection_position = f"{target_position[0]}_C:{target_position[1]}"

        if pre_wash:
            self.transfer_compound(
                specification="liquid",
                compound=pre_wash,
                target_zone=waste_position,
                quantity=pre_wash_volume,
                rinse_volume=0.0
            )

        self.transfer_compound(
            "liquid",
            source_zone=source_zone,
            target_zone=collection_position,
            quantity=volume
        )

        if post_wash:
            self.transfer_compound(
                "liquid",
                compound=post_wash,
                target_zone=collection_position,
                quantity=post_wash_volume,
                rinse_volume=0.0
            )

        self.controller.wait(60)

    def inject_to_hplc(self, source_zone: Union[str, None] = None, filtration: Union[str, None] = None, **kwargs) -> None:
        """
        Performs the injection of a sample to the HPLC-MS by
            - waiting for the instrument status to be "ready"
            - performing the inject and wash operation

        Args:
            source_zone: Source zone to draw solution from.
            filtration: Target zone on the SPE rack given as SPE:$NUMBER.
        """
        if filtration:
            # Quick ugly fix for SPE:$NUMBER -> SPE_D:$NUMBER
            target_position: list = filtration.split(":")
            source_zone = f"{target_position[0]}_D:{target_position[1]}"
            self._mix_liquid_in_well(source_zone)

        print("Waiting for HPLC Status")
        #self.hplc_interface.wait_for_status("ready")
        self._inject_and_wash(source_zone)


    @_operate_if_idle(target_status=None)
    def unmount_all(
            self
    ) -> None:
        """
        Unmounts all tools from the robotic arm.
        """
        self.controller.unmount_all()

    def shutdown_manager(
            self
    ) -> None:
        """
        Shuts down the Manager.app in the AutoSuite Executor.
        """
        self.controller.stop_manager()


    ##################
    # HELPER METHODS #
    ##################

    @_operate_if_idle(target_status="stirring")
    def set_reflux(
            self,
            temperature: float
    ) -> None:
        """
        Starts reflux and stirring at given temperature.

        Args:
            temperature (float)
        """
        self.controller.set_isynth_temperature(state="on", temperature=temperature)
        self.controller.set_isynth_reflux(state="on", temperature=0)
        self.controller.set_stir(stir_zone="ISYNTH", state="on", rpm=300)

    @_operate_if_idle(target_status="idle")
    def stop_reflux(
            self
    ) -> None:
        """
        Stops reflux and stirring.
        """
        self.controller.set_isynth_temperature(state="off")
        self.controller.set_isynth_reflux(state="off")
        self.controller.set_stir(stir_zone="ISYNTH", state="off")

    @_operate_if_idle(target_status=None)
    def _inject_and_wash(
            self,
            source_zone: str
    ):
        """
        Performs the actual injection of a sample to the HPLC injection port.
        Waits for 10 s, then washes the sample loop with a system solvent.

        Args:
            source_zone: Zone name from where to inject to the Injection port.
        """
        self.controller.inject_liquid(source=source_zone, destination="INJECT_L:1", volume=0.5, dst_flow=0.10)
        self.controller.wait(10)
        self._wash_hplc_port()


    @_operate_if_idle(target_status=None)
    def _mix_liquid_in_well(
            self,
            vial: str,
            volume: float = 1.0,
            repeats: int = 3
    ):
        """
        Mixes liquids in a vial by drawing liquid from the vial and dispensing liquid into the same vial afterwards.

        Args:
            vial: Zone definition of the vial.
            volume: Volume to draw and dispense
            repeats: Number of mixing cycles
        """
        # ATTN: Doesn't wash the needle at the current stage.

        for _ in range(repeats):
            self.controller.transfer_liquid(
                source=vial,
                destination=vial,
                volume=volume,
                src_flow=20,
                dst_flow=40,
                dst_td=90,
                rinse_volume=0,
                airgap=0
            )

    def _wash_hplc_port(self, repeats: int = 3):
        """
        I wanna start by apologizing that this method actually exists. The fact that it exists is due to an
        extremely stupid design in AutoSuite. The injection port is washed by transferring liquid to and from a dummy
        well that happens to exist in the exact same physical location as the injection port. After that, the injection
        port is flushed by transferring liquid to and from the injection port itself, followed by a rinse with a system
        solvent.

        Args:
            repeats: Number of wash cycles.
        """
        for _ in range(repeats):
            self.controller.transfer_liquid_bu(
                source="VALVEB:1",
                destination="INJECT_WASH:1",
                volume=0.4,  # That is more than the actual volume of the port, but makes sure to flush everything
                src_flow=20,
                dst_flow=20,
                dst_bu=0,
                rinse_volume=0
            )
            self.controller.transfer_liquid(
                source="INJECT_WASH:1",
                destination="WASTE1:1",
                volume=1.0,  # The transfer volume is larger to make sure that the entire volume is drawn
                src_flow=20,
                dst_flow=40,
                src_bu=0,
                rinse_volume=0
            )

        # flushes the injection port in the 'inject' position
        self.controller.inject_liquid(source="VALVEB:1", destination="INJECT_I:1", volume=0.5, dst_flow=0.10)

        # empties the injection port and washes the needle
        self.controller.transfer_liquid(
            source="INJECT_WASH:1",
            destination="WASTE1:1",
            volume=1.0,  # The transfer volume is larger to make sure that the entire volume is drawn
            src_flow=20,
            dst_flow=40,
            src_bu=0,
            rinse_volume=5
        )