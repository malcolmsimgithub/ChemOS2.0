#!/usr/bin/env python
import numpy as np

class Unit:
    def __init__(self, base_unit_vector, value, system):
        self.vec = base_unit_vector
        self.value = value
        self.system = system

    def __repr__(self):
        return str(self.value)  + self.system.tostr(self.vec)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):
        if isinstance(other, Unit):
            return Unit(self.vec + other.vec, self.value*other.value, self.system)
        else:
            return Unit(self.vec, self.value*other, self.system)

    def __div__(self, other):
        return self.__mul__(other**-1)

    def __add__(self, other):
        uother = self.system.vars['__one'] * other
        if np.all(uother.vec == self.vec):
            return Unit(self.vec, self.value + uother.value, self.system)
        else:
            raise TypeError('Incompatible units')
        
    def __sub__(self, other):
        return self.__add__(other * (-1.0))

    def __pow__(self, other):
        if isinstance(other, float) or isinstance(other, int):
            return Unit(self.vec * other, self.value**other, self.system)
        else:
            raise NotImplementedError('pow not implemented for type '
                                          + str(type(other)))
    def __truediv__(self, other):
        return self * other**(-1)

    def split(self):
        return self.value, Unit(self.vec, 1.0, self.system)

    def express_as(self, expr):
        evexpr = self.system(expr)
        out_val = self/evexpr
        if np.all(out_val.vec==0.0):
            return out_val.value, evexpr
        else:
            return out_val, evexpr

class Units:
    def __init__(self, base_units):
        self.vars = {}
        self.base_units = base_units[:]
        
        for i in range(len(base_units)):
            vec = np.zeros(len(base_units))
            vec[i] = 1.0
            self.vars[base_units[i]] = Unit(vec.copy(), 1.0, self)

        # Make a unitless 1
        self.vars['__one'] = Unit(np.zeros(len(self.base_units)), 1.0, self)

    def tostr(self, vec):
        out = ''
        for i in range(len(self.base_units)):
            if vec[i] == 0.0:
                pass
            else:
                out += ' ' + self.base_units[i] + '^' + str(vec[i])
        return out

    def check(self, name):
        if name in self.vars:
            raise TypeError(name + ' already in unit system.')

    def add_var(self, name, val):
        self.check(name)
        self.vars[name] = self.__call__(val)

    def evaluate_unit_expr(self, expr):
        return eval('__one *' + expr, {}, self.vars)

    def __call__(self, val, unit=None):
        if isinstance(val, float) or isinstance(val, int):
            return self.vars['__one'] * float(val)
        elif isinstance(val, str):
            return self.evaluate_unit_expr(val)
        else:
            raise NotImplementedError()        

def SI(prefixes=True, QM=True, energy=True, em=True, chem=True, extra_base_units=[]):
    s = Units(['m', 's', 'kg', 'A', 'K', 'cd', 'mole'] + extra_base_units)
    s.add_var('pi', np.pi)
    
    if prefixes:
        add_SI_prefix(s)
    if energy:
        add_SI_derived_energy(s)
    if QM:
        add_SI_derived_quantum(s)
    if em:
        add_SI_derived_em(s)
    if chem:
        add_SI_derived_chem(s)
    
    return s

def add_SI_prefix(s):
    s.add_var('peta', 1e15)
    s.add_var('tera', 1e12)
    s.add_var('giga', 1e9)
    s.add_var('mega', 1e6)
    s.add_var('kilo', 1e3)
    s.add_var('deci', 1e-1)
    s.add_var('centi', 1e-2)
    s.add_var('milli', 1e-3)
    s.add_var('micro', 1e-6)
    s.add_var('nano', 1e-9)
    s.add_var('pico', 1e-12)
    s.add_var('femto', 1e-15)
    s.add_var('atto', 1e-18)

def add_SI_derived_energy(s):
    s.add_var('J', 'kg * m**2 / s**2')
    s.add_var('W', 'J/s')
    s.add_var('N', 'J/m')

def add_SI_derived_quantum(s):
    s.add_var('hartree',"4.359744722207185e-18 * J")
    s.add_var('amu', "1.6605390666050e-27 * kg")
    s.add_var('eV', '1.602176620898e-19 * J')
    s.add_var('hbar', '0.658211951440 * eV * femto * s')
    s.add_var('h', 'hbar * 2 * pi')
    s.add_var('ang', '1e-10 * m')
    s.add_var('bohr', "5.2917721090380e-11*m")
    s.add_var('Hz', '1/s')
    s.add_var('charge_e', 'eV * A * s/J')
    s.add_var('mass_e', '9.109383701528e-31*kg')

def add_SI_derived_em(s):
    s.add_var('c', '299792458.0 * m/s')
    s.add_var('C', 'A * s')
    s.add_var('V', 'N*m/C')
    s.add_var('H', 'J/A**2')
    s.add_var('debye', '(1/299792458.0)*10**(-21) * C * m')
    s.add_var('vac_mu', '4 * pi * 1e-7 * H/m')
    s.add_var('vac_e', '1.0/(vac_mu * c**2)')


def add_SI_derived_chem(s):
    s.add_var('litre', '(deci * m)**3')
    s.add_var('molar', 'mole / litre')
    s.add_var('Na', '6.02214076 * 10**23')


def cgs_gaussian(prefixes=True, quantum=True, chem=True, extra_base_units=[]):
    s = Units(['cm', 's', 'g', 'mole'] + extra_base_units)
    s.add_var('pi', np.pi)

    # Various conversion factors
    s.add_var('m', '100 * cm')

    # Energy and forces and such
    s.add_var('dyne', 'g * cm / s**2')
    s.add_var('erg', 'g * cm**2 / s**2')

    # Various units of importance in EM
    s.add_var('statC', 'g**0.5 * cm**1.5 / s')
    s.add_var('c', '299792458.0 * m/s')

    
    if prefixes:
        add_SI_prefix(s)
    if quantum:
        add_cgs_derived_quantum(s)
    if chem:
        add_cgs_derived_chem(s)

    return s

def add_cgs_derived_quantum(s):
    s.add_var('debye', '1e-18 * statC * cm')
    s.add_var('ang', '1e-10 * m')
    s.add_var('h', '6.62607015e-34 * 1e7 * erg * s')
    s.add_var('hbar', 'h/(pi * 2.0)')
    s.add_var('hartree',"4.359744722207185e-18 * 1e7 * erg")

def add_cgs_derived_chem(s):
    s.add_var('litre', '(deci * m)**3')
    s.add_var('molar', 'mole / litre')
    s.add_var('Na', '6.02214076 * 10**23 ')
    s.add_var('bohr', "5.2917721090380e-11*m")
    # s.add_var('eV', '1.602176620898e-19 * J')
    # s.add_var('amu', "1.6605390666050e-27 * kg")
    # s.add_var('Hz', '1/s')
