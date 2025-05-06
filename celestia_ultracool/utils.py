"""General-use functions."""
import numpy as np
import pandas as pd
import re


__all__ = ['app_to_abs_mag', 'calculate_teff', 'parse_name', 'strip_name',
           'get_distance']


def app_to_abs_mag(app_mag: float, dist: float) -> float:
    """Calculate absolute magnitude from apparent magnitude and distance."""
    return app_mag - 5*np.log10(dist/10)


def calculate_teff(lum: float, radius: float) -> tuple[float, str]:
    """Calculate temperature from luminosity and radius, using the Stefan-Boltzmann law."""
    return (10**lum) ** (1/4) * (1/radius) ** (1/2) * 5772, 'from luminosity and radius'


def parse_name(name: str, has_companion: bool=False) -> str:
    """Processes a designation to match the formats used by Celestia."""
    parsed_name = name

    # Remove SIMBAD-specific prefixes
    parsed_name = parsed_name.removeprefix('NAME').removeprefix('V*').removeprefix('EM*')
    parsed_name = re.sub(r'^\*++', '', parsed_name)

    # Trim extra spaces
    parsed_name = parsed_name.strip()
    parsed_name = re.sub(r'\s++', ' ', parsed_name)

    # Separate component letters from the rest of the designation
    if has_companion:
        parsed_name = re.sub(r'(?<=\d)[b-zA-D][a-zB-D]*+$', r' \g<0>', parsed_name)

    # Abbreviate and capitalize Greek letters in Bayer designations
    parsed_name = re.sub(r'^alf', 'ALF', parsed_name)
    parsed_name = re.sub(r'^beta?', 'BET', parsed_name)
    parsed_name = re.sub(r'^chi', 'CHI', parsed_name)
    parsed_name = re.sub(r'^eps', 'EPS', parsed_name)
    parsed_name = re.sub(r'^eta', 'ETA', parsed_name)
    parsed_name = re.sub(r'^kap(?:pa)?', 'KAP', parsed_name)
    parsed_name = re.sub(r'^mu\.02', 'MU2', parsed_name)
    parsed_name = re.sub(r'^pi\.01', 'PI1', parsed_name)
    parsed_name = re.sub(r'^zeta?', 'ZET', parsed_name)

    # For objects with designations from the original Gliese catalog (1-915), write catalog name
    # in full
    parsed_name = re.sub(r'Gl ', 'Gliese ', parsed_name)
    parsed_name = re.sub(r'GJ(?! \d{4})', 'Gliese', parsed_name)

    # Miscellaneous tweaks
    parsed_name = re.sub(r'(?<= \d{3}-)0++', '', parsed_name)
    parsed_name = re.sub(r'^[BC]D ', 'BD', parsed_name)

    return parsed_name


def strip_name(name: str) -> str:
    """Removes the component letter from a designation."""
    return re.sub(r'(?<=\w)\s*+[b-zA-D][a-zB-D]*+$', '', name)


def get_distance(row: pd.Series) -> tuple[float, float, str]:
    """Get non-rounded distance and error from dataset-specific columns."""
    if row.dist_formula_source == 'dist_plx_formula':
        if row.ref_plx_formula in ['Gaia21a', 'Gaia23']:
            plx = row.plx_Gaia
            plx_error = row.plxerr_Gaia
        elif row.ref_plx_formula == 'Best20':
            plx = row.plx_UKIRT
            plx_error = row.plxerr_UKIRT
        elif row.ref_plx_formula == 'Magn20' or pd.isna(row.plx_lit):
            plx = row.plx_formula
            plx_error = row.plxerr_formula
        else:
            plx = row.plx_lit
            plx_error = row.plxerr_lit

        dist = 1000 / plx
        dist_error = plx_error / plx * dist
        dist_note = 'from parallax'
    elif not row.dist_formula_source.startswith('null'):
        dist_column = str(row.dist_formula_source).removesuffix('_young')

        dist = getattr(row, dist_column)
        dist_error = getattr(row, dist_column.replace('dist', 'disterr'))
        dist_note = row.dist_formula_source
    else:
        dist = None
        dist_error = None
        dist_note = 'missing ({})'.format(row.dist_formula_source)

    return dist, dist_error, dist_note
