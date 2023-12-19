import time

# Psuedoparent
class AdvancedPump(object):

    ## must have methods to inherit AdvancedPump
    syringe_volume = None
    def set_valve(self): return
    def get_position(self): return
    def set_position(self): return
    def get_max_steps(self): return

    # AdvancedPump methods
    def convert_port(self, port):
        if isinstance(port, str):
            port = self.ports.get(port)
        if port is not None:
            return port
        else:
            raise ValueError('Wrong ports')

    def draw(self, volume, valve = None): 
        if valve is not None:
            valve = self.convert_port(valve)
            self.set_valve(valve)
        if volume > 0:
            position = self.get_position() + int( (volume / self.syringe_volume) * self.get_max_steps() )
            self.set_position(position)

    def draw_full(self, valve = None): 
        if valve is not None:
            valve = self.convert_port(valve)
            self.set_valve(valve)
        self.set_position(self.get_max_steps())

    def dispense(self, volume, valve = None): 
        if valve is not None:
            valve = self.convert_port(valve)
            self.set_valve(valve)
        if volume > 0:
            position = self.get_position() - int( (volume / self.syringe_volume) * self.get_max_steps() )
            self.set_position(position)

    def dispense_all(self, valve = None):
        if valve is not None:
            valve = self.convert_port(valve)
            self.set_valve(valve)
        self.set_position(0)

    def draw_and_dispense(self, draw_valve, dispense_valve, volume, wait = None, velocity=None):
        "The velocity only temporarily changes here"

        draw_valve     = self.convert_port(draw_valve)
        dispense_valve = self.convert_port(dispense_valve)
        if velocity:
            old_velocity = self.get_velocity()
            self.set_velocity(velocity)

        while volume >= self.syringe_volume:
            volume -= self.syringe_volume
            self.draw_full(draw_valve)
            if wait:
                time.sleep(wait)
            self.dispense_all(dispense_valve)

        self.draw(valve = draw_valve, volume = volume)
        if wait:
            time.sleep(wait)
        self.dispense_all(dispense_valve)

        if velocity:
            self.set_velocity(old_velocity)

    def tmp_velocity(self, velocity):
        class TmpVelocity(object):
            def __enter__(obj):
                obj._old_velocity = self.get_velocity()
                self.set_velocity(velocity)
            def __exit__(obj, *args):
                self.set_velocity(obj._old_velocity)
        return TmpVelocity()

    def __getattr__(self, attr):
        if attr in self.info:
            return self.info[attr]
        else:
            raise ValueError(f'{attr} no exist in obj or obj.info')

