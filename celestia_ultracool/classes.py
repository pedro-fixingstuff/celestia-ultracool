"""Classes used for building single and multiple ultra-cool dwarfs."""

import copy
import numpy as np
from typing import TextIO

from . import consts, bhac15, cond03


class System:
    """A system composed of one or two ultra-cool dwarfs, which can be combined with other
    subsystems to create higher-order hierarchies."""
    def __init__(self, names: list[str], ra: float, dec: float, dist: float, spt_num: float,
                 age: float, infourl: str=None, age_category: str=None, dist_note: str=None,
                 spt_note: str=None) -> None:
        self.names = copy.deepcopy(names)
        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.spt_num = abs(spt_num)
        self.infourl = infourl
        self.age = age
        self.age_category = age_category
        self.dist_note = dist_note
        self.spt_note = spt_note

        self.components = []

        # The evolutionary model grids only go up to 10 Gyr, so clamp the object's age if it exceeds that
        self.log_age = min(np.log10(self.age), 1.0)

    def write(self, stream: TextIO, is_component: bool=False, is_subdwarf: bool=False,
              coord_decimal_digits: int=7, offset_coords: bool=False) -> None:
        """Writes a Celestia STC definition for an ultra-cool dwarf binary/multiple system to a file
        stream."""

        stream.write('\nBarycenter "{}"'.format(':'.join(self.names)))
        stream.write('\n{')
        stream.write('\n\tRA {}'.format(round(self.ra, coord_decimal_digits)))
        stream.write('\n\tDec {}'.format(round(self.dec, coord_decimal_digits)))
        if not is_component and offset_coords:
            stream.write(' # offset from primary coordinates')
        stream.write('\n\tDistance {:.5g}'.format(self.dist))
        if not is_component and self.dist_note is not None:
            stream.write(' # ' + self.dist_note)
        stream.write('\n}\n')

        for component in self.components:
            component.write(stream, True, is_subdwarf, coord_decimal_digits)


class Dwarf:
    """An ultra-cool dwarf."""
    def __init__(self, parent: System=None) -> None:
        self.parent = parent

        self.names = []
        self.ra = parent.ra
        self.dec = parent.dec
        self.spt_num = parent.spt_num
        self.spt_note = parent.spt_note

        self.teff = None
        self.absmag = None
        self.radius = None
        self.teff_note = None
        self.absmag_note = None
        self.radius_note = None
        self.mass = None

    def estimate_teff(self, h_mag: float=0.0, subclass: str='d') -> float:
        """Estimate temperature from empirical relations."""

        # SpT to Teff relation for subdwarfs from Zhang et al. (2018), 2018MNRAS.479.1383Z
        if subclass != 'd' and self.spt_num <= 17:
            self.teff = 3701 - 108.6*self.spt_num + 1.855*self.spt_num**2 - 0.1661*self.spt_num**3

            # The temperature relation applies to [Fe/H] <~ 1.0, that is, esd and usd subclasses.
            # For the sd subclass, an average of the subdwarf and field SpT-Teff relations is taken
            if subclass == 'sd':
                field_teff = self.estimate_teff()
                self.teff = (self.teff + field_teff) / 2
                self.teff_note = 'from spectral type (sd)'
            else:
                self.teff_note = 'from spectral type (esd/usd)'

        elif (self.parent.age <= 0.01 and 1 < self.spt_num < 6) or (self.parent.age < 0.01 and 1 < self.spt_num < 10):
            self.teff = float(np.interp(self.spt_num, consts.SPT_TEFF_YOUNG[0], consts.SPT_TEFF_YOUNG[1]))
            self.teff_note = 'from spectral type (young, <~10 Myr)'

        # M_H to Teff relation (YNG) from Filippazzo et al. (2015), 2015ApJ...810..158F
        elif self.parent.age <= 0.3 and 8.6 < h_mag < 14.5:
            self.teff = -59770 + 23480*h_mag - 3215*h_mag**2 + 190.4*h_mag**3 - 4.167*h_mag**4
            self.teff_note = 'from H-band magnitude (young)'

        # SpT to Teff relation for young objects from Sanghi et al. (2023), 2023ApJ...959...63S
        elif self.parent.age <= 0.3 and 6.0 <= self.spt_num <= 28.0:
            self.teff = (4563.13932 - 228.258451*self.spt_num - 11.1940923*self.spt_num**2
                         + 1.17847855*self.spt_num**3 - 0.0238265571*self.spt_num**4)
            self.teff_note = 'from spectral type (young)'

        # M_H to Teff relation from Kirkpatrick et al. (2021), 2021ApJS..253....7K
        elif 9.5 <= h_mag <= 25.0:
            self.teff = 12516 - 1566.6*h_mag + 67.502*h_mag**2 - 0.92430*h_mag**3 - 0.0019530*h_mag**4
            self.teff_note = 'from H-band magnitude (field)'

        # Piecewise SpT to Teff relation from Kirkpatrick et al. (2021)
        # While the paper states that it only is valid from type L0, it can be seen from both
        # Table 14 and Figure 20 that it also applies to late-M dwarfs
        elif 6.0 <= self.spt_num <= 32.0:
            if self.spt_num < 18.75:
                self.teff = 2237.5 - 144.96*(self.spt_num-10) + 4.0301*(self.spt_num-10)**2
            elif 18.75 <= self.spt_num < 24.75:
                self.teff = 1437.9 - 18.309*(self.spt_num-10)
            else:
                self.teff = 5141.3 - 368.65*(self.spt_num-10) + 6.7301*(self.spt_num-10)**2
            self.teff_note = 'from spectral type (field)'

        else:
            self.teff = float(np.interp(self.spt_num, consts.SPT_TEFF[0], consts.SPT_TEFF[1]))
            self.teff_note = 'from spectral type'

        return self.teff

    def estimate_absmag(self, lum: float=None) -> None:
        """Estimate visual absolute magnitude from evolutionary model grids."""
        if lum is not None:
            if lum > bhac15.min_lum(self.parent.log_age):
                self.absmag = float(bhac15.t_l_mv_interp(self.parent.log_age, lum))
                self.absmag_note = 'BHAC'
            else:
                self.absmag = float(cond03.t_l_mv_interp(self.parent.log_age, lum))
                self.absmag_note = 'COND'
        else:
            log_teff = np.log10(self.teff)

            if self.teff > bhac15.min_teff(self.parent.log_age):
                self.absmag = float(bhac15.t_teff_mv_interp(self.parent.log_age, log_teff))
                self.absmag_note = 'BHAC'
            else:
                self.absmag = float(cond03.t_teff_mv_interp(self.parent.log_age, log_teff))
                self.absmag_note = 'COND'

    def estimate_radius(self, lum: float=None) -> None:
        """Estimate radius from evolutionary model grids."""
        if lum is not None:
            if lum > bhac15.min_lum(self.parent.log_age):
                log_radius = bhac15.t_l_r_interp(self.parent.log_age, lum)
                radius_note = 'BHAC'
            else:
                log_radius = cond03.t_l_r_interp(self.parent.log_age, lum)
                radius_note = 'COND'
            self.radius = 10 ** log_radius * consts.SOLAR_RADIUS
            self.radius_note = radius_note
        else:
            log_teff = np.log10(self.teff)

            if self.teff > bhac15.min_teff(self.parent.log_age):
                log_radius = bhac15.t_teff_r_interp(self.parent.log_age, log_teff)
                radius_note = 'BHAC'
            else:
                log_radius = cond03.t_teff_r_interp(self.parent.log_age, log_teff)
                radius_note = 'COND'
            self.radius = 10 ** log_radius * consts.SOLAR_RADIUS
            self.radius_note = radius_note

    def estimate_mass(self, lum: float=None) -> None:
        """Estimate mass from evolutionary model grids."""
        if lum is not None:
            if lum > bhac15.min_lum(self.parent.log_age):
                self.mass = 10 ** bhac15.t_l_m_interp(self.parent.log_age, lum)
            else:
                self.mass = 10 ** cond03.t_l_m_interp(self.parent.log_age, lum)
        else:
            log_teff = np.log10(self.teff)

            if self.teff > bhac15.min_teff(self.parent.log_age):
                self.mass = 10 ** bhac15.t_teff_m_interp(self.parent.log_age, log_teff)
            else:
                self.mass = 10 ** cond03.t_teff_m_interp(self.parent.log_age, log_teff)

    def write(self, stream: TextIO, is_component: bool=False, is_subdwarf: bool=False,
              coord_decimal_digits: int=7, offset_coords: bool=False) -> None:
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
        stream.write('\n\t# Assigned age: {} Gyr'.format(self.parent.age))
        if self.parent.age_category is not None:
            stream.write(' ({})'.format(self.parent.age_category))
        stream.write('\n\tRA {}'.format(round(self.ra, coord_decimal_digits)))
        stream.write('\n\tDec {}'.format(round(self.dec, coord_decimal_digits)))
        if not is_component and offset_coords:
            stream.write(' # offset from primary coordinates')
        stream.write('\n\tDistance {:.5g}'.format(self.parent.dist))
        if not is_component and self.parent.dist_note is not None:
            stream.write(' # ' + self.parent.dist_note)
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
        if self.parent.infourl is not None:
            stream.write('\n\tInfoURL "{}"'.format(self.parent.infourl))
        stream.write('\n}\n')
