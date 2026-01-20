"""Classes used for building single and multiple ultra-cool dwarfs."""

import copy
import numpy as np
from typing import TextIO

from . import consts
from .bhac15 import BhacInterpolator
from .cond03 import CondInterpolator


class EvoInterpolator:
    """Interpolate physical properties from different evolutionary models."""
    def __init__(self) -> None:
        self.bhac_interp = BhacInterpolator()
        self.cond_interp = CondInterpolator()

    def interpolate_absmag(self, log_age: float, lbol: float, teff: float) -> tuple[float, str]:
        """Interpolate visual absolute magnitude from evolutionary model grids."""
        if lbol is not None:
            if lbol > self.bhac_interp.min_lum(log_age):
                return self.bhac_interp.interpolate_t_l_mv(log_age, lbol), 'BHAC'
            else:
                return self.cond_interp.interpolate_t_l_mv(log_age, lbol), 'COND'
        else:
            log_teff = np.log10(teff)

            if teff > self.bhac_interp.min_teff(log_age):
                return self.bhac_interp.interpolate_t_teff_mv(log_age, log_teff), 'BHAC'
            else:
                return self.cond_interp.interpolate_t_teff_mv(log_age, log_teff), 'COND'

    def interpolate_radius(self, log_age: float, lbol: float, teff: float) -> tuple[float, str]:
        """Interpolate radius from evolutionary model grids."""
        if lbol is not None:
            if lbol > self.bhac_interp.min_lum(log_age):
                return self.bhac_interp.interpolate_t_l_r(log_age, lbol), 'BHAC'
            else:
                return self.cond_interp.interpolate_t_l_r(log_age, lbol), 'COND'
        else:
            log_teff = np.log10(teff)

            if teff > self.bhac_interp.min_teff(log_age):
                return self.bhac_interp.interpolate_t_teff_r(log_age, log_teff), 'BHAC'
            else:
                return self.cond_interp.interpolate_t_teff_r(log_age, log_teff), 'COND'

    def interpolate_mass(self, log_age: float, lbol: float, teff: float) -> float:
        """Interpolate mass from evolutionary model grids."""
        if lbol is not None:
            if lbol > self.bhac_interp.min_lum(log_age):
                return self.bhac_interp.interpolate_t_l_m(log_age, lbol)
            else:
                return self.cond_interp.interpolate_t_l_m(log_age, lbol)
        else:
            log_teff = np.log10(teff)

            if teff > self.bhac_interp.min_teff(log_age):
                return self.bhac_interp.interpolate_t_teff_m(log_age, log_teff)
            else:
                return self.cond_interp.interpolate_t_teff_m(log_age, log_teff)


class Orbit:
    """A Keplerian orbit of a binary system."""
    def __init__(self, period: float=None, sma: float=None, ecc: float=0.0, inc: float=0.0,
                 node: float=0.0, arg_peri: float=0.0, mean_anomaly: float=0.0, epoch: float=consts.J2000) -> None:
        self.period = period
        self.sma = sma
        self.sma_1 = sma
        self.sma_2 = sma
        self.ecc = ecc
        self.inc = inc
        self.node = node
        self.arg_peri = arg_peri
        self.mean_anomaly = mean_anomaly
        self.epoch = epoch

    def write(self, component_idx: int, pri_mass: float, sec_mass: float, measured_masses: bool=False) -> str:
        """Writes an orbit definition as Celestia STC parameters."""
        output = ''

        if component_idx == 0:  # primary
            sma = self.sma_1
            arg_peri = (self.arg_peri + 180) % 360
        elif component_idx == 1:  # secondary
            sma = self.sma_2
            arg_peri = self.arg_peri

        output += '\n\tEllipticalOrbit\n\t{'
        output += '\n\t\tPeriod {:g}'.format(self.period)
        output += '\n\t\tSemiMajorAxis {:.4g}'.format(sma)
        if measured_masses:
            output += ' # mass ratio {:.3f}:{:.3f}'.format(pri_mass, sec_mass)
        else:
            output += ' # mass ratio approximate'
        if self.ecc:
            output += '\n\t\tEccentricity {}'.format(self.ecc)
        if self.inc:
            output += '\n\t\tInclination {}'.format(round(self.inc, 2))
        if self.node:
            output += '\n\t\tAscendingNode {}'.format(round(self.node, 2))
        if arg_peri:
            output += '\n\t\tArgOfPericenter {}'.format(round(arg_peri, 2))
        if self.mean_anomaly:
            output += '\n\t\tMeanAnomaly {}'.format(round(self.mean_anomaly, 2))
        if self.epoch != consts.J2000:
            output += '\n\t\tEpoch {}'.format(round(self.epoch, 6))
        output += '\n\t}'

        return output


class UltracoolSystem:
    """A system composed of one or two ultra-cool dwarfs, which can be combined with other
    subsystems to create higher-order hierarchies."""
    def __init__(self, names: list[str], ra: float, dec: float, dist: float, spt_num: float, age: float, *,
                 subclass: str='d', infourl: str=None, age_category: str=None, dist_note: str=None, spt_note: str=None,
                 parent: 'UltracoolSystem'=None, components: list['UltracoolDwarf'] | list['UltracoolSystem']=None,
                 mass: float=None, orbit: Orbit=None) -> None:
        self.names = copy.deepcopy(names)
        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.spt_num = abs(spt_num)
        self.age = age
        self.subclass = subclass
        self.infourl = infourl
        self.age_category = age_category
        self.dist_note = dist_note
        self.spt_note = spt_note

        self.parent = parent
        self.components = [] if components is None else components
        self.mass = mass
        self.orbit = orbit

        # The evolutionary model grids only go up to 10 Gyr, so clamp the object's age if it exceeds that
        self.log_age = min(np.log10(self.age), 1.0)

    def write(self, stream: TextIO, *, is_component: bool=False, is_subdwarf: bool=False,
              coord_decimal_digits: int=7, offset_coords: bool=False, measured_masses: bool=False) -> None:
        """Writes a Celestia STC definition for an ultra-cool dwarf binary/multiple system to a file
        stream."""

        stream.write('\nBarycenter "{}"'.format(':'.join(self.names)))
        stream.write('\n{')

        if self.parent is not None and self.parent.orbit is not None:
            stream.write('\n\tOrbitBarycenter "{}"'.format(self.parent.names[0]))

            index = self.parent.components.index(self)
            if index == 0:  # primary
                pri_mass = self.mass
                sec_mass = self.parent.components[1].mass
            else:  # secondary
                pri_mass = self.parent.components[0].mass
                sec_mass = self.mass

            stream.write(self.parent.orbit.write(index, pri_mass, sec_mass, measured_masses))
        else:
            stream.write('\n\tRA {}'.format(round(self.ra, coord_decimal_digits)))
            stream.write('\n\tDec {}'.format(round(self.dec, coord_decimal_digits)))
            if not is_component and offset_coords:
                stream.write(' # offset from primary coordinates')
            stream.write('\n\tDistance {:.5g}'.format(self.dist))
            if not is_component and self.dist_note is not None:
                stream.write(' # ' + self.dist_note)

        stream.write('\n}\n')

        for component in self.components:
            component.write(stream, is_component=True, is_subdwarf=is_subdwarf,
                            coord_decimal_digits=coord_decimal_digits, measured_masses=measured_masses)


class UltracoolDwarf:
    """An ultra-cool dwarf."""
    def __init__(self, parent: UltracoolSystem=None, mags: dict[str, float]=None,
                 lbol: float=None, teff: float=None, absmag: float=None, radius: float=None, mass: float=None,
                 teff_note: str=None, absmag_note: str=None, radius_note: str=None, properties_note: str=None) -> None:
        self.parent = parent

        self.names = []
        self.ra = parent.ra
        self.dec = parent.dec
        self.spt_num = parent.spt_num
        self.spt_note = parent.spt_note

        self.mags = {} if mags is None else mags

        self.lbol = lbol
        self.teff = teff
        self.absmag = absmag
        self.radius = radius
        self.mass = mass
        self.teff_note = teff_note
        self.absmag_note = absmag_note
        self.radius_note = radius_note
        self.properties_note = properties_note

    def estimate_lbol_teff(self) -> None:
        """Estimate luminosity or temperature from empirical relations."""

        # SpT to Teff relation for subdwarfs from Zhang et al. (2018), 2018MNRAS.479.1383Z
        if self.parent.subclass != 'd' and self.spt_num <= 17:
            self.teff = 3701 - 108.6*self.spt_num + 1.855*self.spt_num**2 - 0.1661*self.spt_num**3

            # The temperature relation applies to [Fe/H] <~ 1.0, that is, esd and usd subclasses.
            # For the sd subclass, an average of the subdwarf and field SpT-Teff relations is taken
            if self.parent.subclass == 'sd':
                field_teff = float(np.interp(self.spt_num, consts.SPT_TEFF[0], consts.SPT_TEFF[1]))
                self.teff = (self.teff + field_teff) / 2
                self.properties_note = 'estimated from spectral type (sd)'
            else:
                self.properties_note = 'estimated from spectral type (esd/usd)'

        # Relations for young objects
        elif self.parent.age <= 0.3:
            if (self.parent.age <= 0.01 and self.spt_num <= 9):
                self.teff = float(np.interp(self.spt_num, consts.SPT_TEFF_YOUNG[0], consts.SPT_TEFF_YOUNG[1]))
                self.properties_note = 'estimated from spectral type (young, <~10 Myr)'

            # Absolute magnitude to Lbol polynomials from Sanghi et al. (2023), 2023ApJ...959...63S
            elif 7.7 <= self.mags.get('H_MKO', 0) <= 17.0:  # rms = 0.057 dex
                self.lbol = 7.25956089e-1 - 3.88884874e-1*self.mags['H_MKO']
                self.properties_note = 'estimated from MKO H-band (young)'
            elif 7.3 <= self.mags.get('K_MKO', 0) <= 16.7:  # rms = 0.057 dex
                self.lbol = (-2.06743904e1 + 7.40114284*self.mags['K_MKO'] - 1.03611814*self.mags['K_MKO']**2
                             + 5.90347236e-2*self.mags['K_MKO']**3 - 1.21640941e-3*self.mags['K_MKO']**4)
                self.properties_note = 'estimated from MKO K-band (young)'
            elif 7.3 <= self.mags.get('Ks', 0) <= 16.8:  # rms = 0.057 dex
                self.lbol = (-1.68245519e1 + 5.94124964*self.mags['Ks'] - 8.33426219e-1*self.mags['Ks']**2
                             + 4.68755666e-2*self.mags['Ks']**3 - 9.51087115e-4*self.mags['Ks']**4)
                self.properties_note = 'estimated from 2MASS Ks-band (young)'
            elif 7.7 <= self.mags.get('H', 0) <= 17.0:  # rms = 0.072 dex
                self.lbol = 6.58828628e-1 - 3.84611365e-1*self.mags['H']
                self.properties_note = 'estimated from 2MASS H-band (young)'
            elif 6.6 <= self.mags.get('W1', 0) <= 16.1:  # rms = 0.116 dex
                self.lbol = (-3.70291470e+1 + 1.38286472e+1*self.mags['W1'] - 1.95883111*self.mags['W1']**2
                             + 1.15867819e-1*self.mags['W1']**3 - 2.48769429e-3*self.mags['W1']**4)
                self.properties_note = 'estimated from WISE W1-band (young)'
            elif 11.0 <= self.mags.get('z', 0) <= 19.3:  # rms = 0.121 dex
                self.lbol = 6.92111297e-1 - 2.88226403e-1*self.mags['z']
                self.properties_note = 'estimated from PS1 z-band (young)'
            elif 11.0 <= self.mags.get('W2', 0) <= 19.3:  # rms = 0.158 dex
                self.lbol = 1.05360011 - 4.84834949e-1*self.mags['W2']
                self.properties_note = 'estimated from WISE W2-band (young)'
            # ML relations
            elif self.spt_num < 20 and 8.2 <= self.mags.get('J_MKO', 0) <= 15.8:  # rms = 0.080 dex
                self.lbol = 4.06864509e-1 - 3.36361777e-1*self.mags['J_MKO']
                self.properties_note = 'estimated from MKO J-band (young)'
            elif self.spt_num < 20 and 8.2 <= self.mags.get('J', 0) <= 15.8:  # rms = 0.081 dex
                self.lbol = 4.08756088e-1 - 3.34689569e-1*self.mags['J']
                self.properties_note = 'estimated from 2MASS J-band (young)'
            elif self.spt_num < 20 and 10.3 <= self.mags.get('y', 0) <= 17.8:  # rms = 0.104 dex
                self.lbol = 3.28920351e-1 - 2.80644117e-1*self.mags['y']
                self.properties_note = 'estimated from PS1 y-band (young)'
            # T relations
            elif 20 <= self.spt_num < 30 and 16.5 <= self.mags.get('y', 0) <= 19.5:  # rms = 0.075 dex
                self.lbol = 6.92111297e-1 - 2.88226403e-1*self.mags['y']
                self.properties_note = 'estimated from PS1 y-band (young)'
            elif 20 <= self.spt_num < 30 and 14.0 <= self.mags.get('J_MKO', 0) <= 16.7:  # rms = 0.097 dex
                self.lbol = 8.31005313e-1 - 3.94047446e-1*self.mags['J_MKO']
                self.properties_note = 'estimated from MKO J-band (young)'
            elif 20 <= self.spt_num < 30 and 14.1 <= self.mags.get('J', 0) <= 16.7:  # rms = 0.108 dex
                self.lbol = 9.53221827e-1 - 3.95259421e-1*self.mags['J']
                self.properties_note = 'estimated from 2MASS J-band (young)'

            # SpT to Teff relation for young objects from Sanghi et al. (2023)
            elif 6.0 <= self.spt_num <= 28.0:
                self.teff = (4.56313932e3 - 2.28258451e2*self.spt_num - 1.11940923e1*self.spt_num**2
                             + 1.17847855*self.spt_num**3 - 2.38265571e-2*self.spt_num**4)
                self.properties_note = 'estimated from spectral type (young)'

        # Field relations, also used for young objects if they span ranges not covered above
        if self.lbol is None and self.teff is None:
            # M_W2 to Teff relation from Leggett et al. (2025), 2025ApJ...991..193L
            if 12.8 <= self.mags.get('W2', 0) <= 17.5:
                self.teff = (258091 - 64036.1*self.mags['W2'] + 5979.92*self.mags['W2']**2
                             - 248.49*self.mags['W2']**3 + 3.87238*self.mags['W2']**4)
                self.properties_note = 'estimated from WISE W2-band'

            # Absolute magnitude to Lbol polynomials from Sanghi et al. (2023)
            elif 8.5 <= self.mags.get('H_MKO', 0) <= 18.8:  # rms = 0.053 dex
                self.lbol = 6.17668495e-1 - 3.80907000e-1*self.mags['H_MKO']
                self.properties_note = 'estimated from MKO H-band'
            elif 8.1 <= self.mags.get('K_MKO', 0) <= 19.0:  # rms = 0.064 dex
                self.lbol = (-8.68942422 + 3.17477734*self.mags['K_MKO'] - 4.83389959e-1*self.mags['K_MKO']**2
                             + 2.73739037e-2*self.mags['K_MKO']**3 - 5.47040525e-4*self.mags['K_MKO']**4)
                self.properties_note = 'estimated from MKO K-band'
            elif 8.5 <= self.mags.get('H', 0) <= 18.1:  # rms = 0.066 dex
                self.lbol = 6.93323688e-1 - 3.89636901e-1*self.mags['H']
                self.properties_note = 'estimated from 2MASS H-band'
            elif 11.0 <= self.mags.get('z', 0) <= 22.7:  # rms = 0.083 dex
                self.lbol = 8.78235147e-1 - 3.06670163e-1*self.mags['z']
                self.properties_note = 'estimated from PS1 z-band'
            elif 7.3 <= self.mags.get('Ks', 0) <= 18.5:  # rms = 0.083 dex
                self.lbol = (-2.50181488e1 + 8.16617507*self.mags['Ks'] - 1.04402978*self.mags['Ks']**2
                             + 5.48196484e-2*self.mags['Ks']**3 - 1.04209064e-3*self.mags['Ks']**4)
                self.properties_note = 'estimated from 2MASS Ks-band'
            elif 7.6 <= self.mags.get('W1', 0) <= 17.1:  # rms = 0.103 dex
                self.lbol = (-5.22934625e+1 + 1.85121627e+1*self.mags['W1'] - 2.47203555*self.mags['W1']**2
                             + 1.39599009e-1*self.mags['W1']**3 - 2.87497133e-3*self.mags['W1']**4)
                self.properties_note = 'estimated from WISE W1-band'
            elif 7.3 <= self.mags.get('W2', 0) <= 14.1:  # rms = 0.125 dex
                self.lbol = 1.75395276 - 5.42818335e-1*self.mags['W2']
                self.properties_note = 'estimated from WISE W2-band'
            # ML relations
            elif self.spt_num < 20 and 9.0 <= self.mags.get('J_MKO', 0) <= 16.0:  # rms = 0.055 dex
                self.lbol = 6.00343959e-1 - 3.57204309e-1*self.mags['J_MKO']
                self.properties_note = 'estimated from MKO J-band'
            elif self.spt_num < 20 and 9.1 <= self.mags.get('J', 0) <= 16.1:  # rms = 0.057 dex
                self.lbol = 5.53301206e-1 - 3.51344336e-1*self.mags['J']
                self.properties_note = 'estimated from 2MASS J-band'
            elif self.spt_num < 20 and 10.5 <= self.mags.get('y', 0) <= 18.5:  # rms = 0.074 dex
                self.lbol = 5.55740568e-1 - 3.03933139e-1*self.mags['y']
                self.properties_note = 'estimated from PS1 y-band'
            # T relations
            elif 20 <= self.spt_num < 30 and 15.0 <= self.mags.get('y', 0) <= 21.0:  # rms = 0.066 dex
                self.lbol = 1.37035364 - 3.66205917e-1*self.mags['y']
                self.properties_note = 'estimated from PS1 y-band'
            elif 20 <= self.spt_num < 30 and 12.8 <= self.mags.get('J', 0) <= 18.4:  # rms = 0.101 dex
                self.lbol = 8.80328596e-1 - 3.93121219e-1*self.mags['J']
                self.properties_note = 'estimated from 2MASS J-band'
            elif 20 <= self.spt_num < 30 and 12.6 <= self.mags.get('J_MKO', 0) <= 18.4:  # rms = 0.108 dex
                self.lbol = 9.54537740e-1 - 4.04863108e-1*self.mags['J_MKO']
                self.properties_note = 'estimated from MKO J-band'

            # SpT to Lbol relation
            elif 6.0 <= self.spt_num <= 29.0:
                self.lbol = (2.38575339 - 2.34765187*self.spt_num + 4.10403429e-1*self.spt_num**2
                             - 3.75870300e-2*self.spt_num**3 + 1.79418353e-3*self.spt_num**4
                             - 4.22643837e-5*self.spt_num**5 + 3.85072198e-7*self.spt_num**6)
                self.properties_note = 'estimated from spectral type'

            else:
                self.teff = float(np.interp(self.spt_num, consts.SPT_TEFF[0], consts.SPT_TEFF[1]))
                self.properties_note = 'estimated from spectral type'

    def estimate_properties(self, interp: EvoInterpolator):
        """Estimate fundamental properties from empirical relations and evolutionary models."""
        if self.lbol is None:
            self.estimate_lbol_teff()
        self.absmag, self.absmag_note = interp.interpolate_absmag(self.parent.log_age, self.lbol, self.teff)
        if self.radius is None:
            self.radius, self.radius_note = interp.interpolate_radius(self.parent.log_age, self.lbol, self.teff)
            self.radius *= consts.SOLAR_RADIUS
            if self.teff is None:
                # Stefan-Boltzmann law
                self.teff = (10**self.lbol) ** (1/4) * (consts.SOLAR_RADIUS/self.radius) ** (1/2) * consts.SOLAR_TEFF
        if self.mass is None:
            self.mass = interp.interpolate_mass(self.parent.log_age, self.lbol, self.teff)

    def write(self, stream: TextIO, *, is_component: bool=False, is_subdwarf: bool=False,
              coord_decimal_digits: int=7, offset_coords: bool=False, measured_masses: bool=False) -> None:
        """Writes a Celestia STC definition for an ultra-cool dwarf to a file stream."""

        # Extract spectral class from SpT number
        if self.spt_num < 10:
            spt = 'M'
        elif self.spt_num < 20:
            spt = 'L'
        elif self.spt_num < 30:
            spt = 'T'
        else:
            spt = 'Y'

        # Append the subtype
        spt += '{:g}'.format(self.spt_num % 10)

        # Negative values indicate subdwarfs
        if is_subdwarf:
            spt = 'sd' + spt

        stream.write('\n"{}"'.format(':'.join(self.names)))
        stream.write('\n{')
        stream.write('\n\t# Assigned age: {:g} Gyr'.format(self.parent.age))
        if self.parent.age_category is not None:
            stream.write(' ({})'.format(self.parent.age_category))
        if self.properties_note is not None:
            stream.write('\n\t# Fundamental properties {}'.format(self.properties_note))

        if self.parent.orbit is None:
            stream.write('\n\tRA {}'.format(round(self.ra, coord_decimal_digits)))
            stream.write('\n\tDec {}'.format(round(self.dec, coord_decimal_digits)))
            if not is_component and offset_coords:
                stream.write(' # offset from primary coordinates')
            stream.write('\n\tDistance {:.5g}'.format(self.parent.dist))
            if not is_component and self.parent.dist_note is not None:
                stream.write(' # ' + self.parent.dist_note)
        else:
            stream.write('\n\tOrbitBarycenter "{}"'.format(self.parent.names[0]))

        stream.write('\n\tSpectralType "{}"'.format(spt))
        if self.spt_note is not None:
            stream.write(' # ' + self.spt_note)
        stream.write('\n\tTemperature {}'.format(round(self.teff)))
        if self.teff_note is not None:
            stream.write(' # ' + self.teff_note)
        stream.write('\n\tAbsMag {}'.format(round(self.absmag, 2)))
        if self.absmag_note is not None:
            stream.write(' # ' + self.absmag_note)
        stream.write('\n\tRadius {}'.format(round(self.radius)))
        if self.radius_note is not None:
            stream.write(' # ' + self.radius_note)

        if self.parent.orbit is not None:
            index = self.parent.components.index(self)
            if index == 0:  # primary
                pri_mass = self.mass
                sec_mass = self.parent.components[1].mass
            else:  # secondary
                pri_mass = self.parent.components[0].mass
                sec_mass = self.mass

            stream.write(self.parent.orbit.write(index, pri_mass, sec_mass, measured_masses))

        if self.parent.infourl is not None:
            stream.write('\n\tInfoURL "{}"'.format(self.parent.infourl))
        stream.write('\n}\n')
