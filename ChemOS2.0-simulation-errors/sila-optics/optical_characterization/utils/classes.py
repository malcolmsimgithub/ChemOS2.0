


class Absorption(object):
    def __init__(self, data = None):
        self.attrs = ['reference', 'sample', 'transmittance', 'absorbance', 'max', 'lambda_max', 'end_index', 'end_wavelength']
        self._data = {attr : None for attr in self.attrs}

        if data:
            for key, val in data.items():
                if key in self._data.keys():
                    self._data['key'] = val

    @property
    def reference(self):
        return self._data['reference']

    @property
    def sample(self):
        return self._data['sample']

    @property
    def transmittance(self):
        return self._data['transmittance']

    @property
    def absorbance(self):
        return self._data['absorbance']

    @property
    def max(self):
        return self._data['max']

    @property
    def lambda_max(self):
        return self._data['lambda_max']

    @property
    def end_index(self):
        return self._data['end_index']

    @property
    def end_wavelength(self):
        return self._data['end_wavelength']

    @reference.setter
    def reference(self, value):
        self._data['reference'] = value

    @sample.setter
    def sample(self, value):
        self._data['sample'] = value

    @transmittance.setter
    def transmittance(self, value):
        self._data['transmittance'] = value

    @absorbance.setter
    def absorbance(self, value):
        self._data['absorbance'] = value

    @max.setter
    def max(self, value):
        self._data['max'] = value

    @lambda_max.setter
    def lambda_max(self, value):
        self._data['lambda_max'] = value

    @end_index.setter
    def end_index(self, value):
        self._data['end_index'] = value

    @end_wavelength.setter
    def end_wavelength(self, value):
        self._data['end_wavelength'] = value




class UV(object):
    def __init__(self, data = None):
        self.attrs = ['reference', 'absorption', 'absorbance_time_trace', 'absorbance_maintenance', 'degradation_rate', 'absorbance_maintenance(at_1min)', 'wavelength']
        self._data = {attr : None for attr in self.attrs}

        if data:
            for key, val in data.items():
                if key in self._data.keys():
                    self._data['key'] = val

    @property
    def reference(self):
        return self._data['reference']

    @reference.setter
    def reference(self, value):
        self._data['reference'] = value

    @property
    def absorption(self):
        return self._data['absorption']

    @absorption.setter
    def absorption(self, value):
        self._data['absorption'] = value

    @property
    def absorbance_time_trace(self):
        return self._data['absorbance_time_trace']

    @absorbance_time_trace.setter
    def absorbance_time_trace(self, value):
        self._data['absorbance_time_trace'] = value

    @property
    def absorbance_maintenance(self):
        return self._data['absorbance_maintenance']

    @absorbance_maintenance.setter
    def absorbance_maintenance(self, value):
        self._data['absorbance_maintenance'] = value

    @property
    def degradation_rate(self):
        return self._data['degradation_rate']

    @degradation_rate.setter
    def degradation_rate(self, value):
        self._data['degradation_rate'] = value

    @property
    def absorbance_maintenance_at_1min(self):
        return self._data['absorbance_maintenance(at_1min)']

    @absorbance_maintenance_at_1min.setter
    def absorbance_maintenance_at_1min(self, value):
        self._data['absorbance_maintenance(at_1min)'] = value

    @property
    def wavelength(self):
        return self._data['wavelength']

    @wavelength.setter
    def wavelength(self, value):
        self._data['wavelength'] = value

class Emission(object):
    def __init__(self, data = None):
        self.attrs = ['energy', 'photons', 'freq_spectrum', 'gain_spectrum', 'relative_QY', 'max', 'lambda_max', 'time_trace', 'int_time_trace', \
                 'maintenance', 'degradation_rate', 'maintenance(at_1min)', 'exposure', 'max_gain_factor', 'max_gain_wavelength']
        self._data = {attr : None for attr in self.attrs}

        if data:
            for key, val in data.items():
                if key in self._data.keys():
                    self._data['key'] = val


    @property
    def energy(self):
        return self._data['energy']

    @energy.setter
    def energy(self, value):
        self._data['energy'] = value

    @property
    def photons(self):
        return self._data['photons']

    @photons.setter
    def photons(self, value):
        self._data['photons'] = value

    @property
    def freq_spectrum(self):
        return self._data['freq_spectrum']

    @freq_spectrum.setter
    def freq_spectrum(self, value):
        self._data['freq_spectrum'] = value

    @property
    def gain_spectrum(self):
        return self._data['gain_spectrum']

    @gain_spectrum.setter
    def gain_spectrum(self, value):
        self._data['gain_spectrum'] = value

    @property
    def relative_QY(self):
        return self._data['relative_QY']

    @relative_QY.setter
    def relative_QY(self, value):
        self._data['relative_QY'] = value

    @property
    def max(self):
        return self._data['max']

    @max.setter
    def max(self, value):
        self._data['max'] = value

    @property
    def lambda_max(self):
        return self._data['lambda_max']

    @lambda_max.setter
    def lambda_max(self, value):
        self._data['lambda_max'] = value

    @property
    def time_trace(self):
        return self._data['time_trace']

    @time_trace.setter
    def time_trace(self, value):
        self._data['time_trace'] = value

    @property
    def int_time_trace(self):
        return self._data['int_time_trace']

    @int_time_trace.setter
    def int_time_trace(self, value):
        self._data['int_time_trace'] = value

    @property
    def maintenance(self):
        return self._data['maintenance']

    @maintenance.setter
    def maintenance(self, value):
        self._data['maintenance'] = value

    @property
    def degradation_rate(self):
        return self._data['degradation_rate']

    @degradation_rate.setter
    def degradation_rate(self, value):
        self._data['degradation_rate'] = value

    @property
    def maintenance_at_1min(self):
        return self._data['maintenance(at_1min)']

    @maintenance_at_1min.setter
    def maintenance_at_1min(self, value):
        self._data['maintenance(at_1min)'] = value

    @property
    def exposure(self):
        return self._data['exposure']

    @exposure.setter
    def exposure(self, value):
        self._data['exposure'] = value

    @property
    def max_gain_factor(self):
        return self._data['max_gain_factor']

    @max_gain_factor.setter
    def max_gain_factor(self, value):
        self._data['max_gain_factor'] = value

    @property
    def max_gain_wavelength(self):
        return self._data['max_gain_wavelength']

    @max_gain_wavelength.setter
    def max_gain_wavelength(self, value):
        self._data['max_gain_wavelength'] = value