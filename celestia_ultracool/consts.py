"""General-use constants."""

# Gaia DR3 epoch, used by Celestia stars.dat
DEFAULT_EPOCH = 2016.0

# Mean values from the main table
PLACEHOLDER_DIST = 180  # ly, rounded
PLACEHOLDER_SPT_NUM = 15  # rounded

# SpT-Teff scale from Pecaut and Mamajek (2013), available at
# https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
SPT_TEFF = ((
        0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 11, 12,
        13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 24.5, 25, 25.5, 26, 27, 27.5, 28, 28.5, 29,
        29.5, 30, 30.5, 31, 31.5, 32, 34), (
        3850, 3770, 3660, 3620, 3560, 3470, 3430, 3270, 3210, 3110, 3060, 2930, 2810, 2740, 2680,
        2630, 2570, 2420, 2380, 2350, 2270, 2160, 2060, 1920, 1870, 1710, 1550, 1530, 1420, 1370,
        1255, 1240, 1220, 1200, 1180, 1170, 1160, 1040, 950, 825, 750, 680, 600, 560, 510, 450, 400,
        360, 325, 320, 250))

# SpT-Teff scale for M-type T Tauri stars from Herczeg and Hillenbrand (2014), 2014ApJ...786...97H
SPT_TEFF_YOUNG = ((0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                  (3900, 3720, 3560, 3410, 3190, 2980, 2860, 2770, 2670, 2570))

J2000 = 2451545.0

# Mean solar and Jovian radii
SOLAR_RADIUS = 695700.0
JUPITER_RADIUS = 69911.0

SOLAR_TEFF = 5772

# Mass conversion factor
MSUN_TO_MJUP = 1.3271244e20 / 1.2668653e17

YEAR_TO_DAY = 365.25
PC_TO_LY = 3.261563777

HEADER = """# Catalog of ultra-cool dwarfs, based on the UltracoolSheet (https://zenodo.org/records/10573247)
# maintained by Will Best, Trent Dupuy, Michael Liu, Aniket Sanghi, Rob Siverd and Zhoujian Zhang.

# Coordinates are normalized to the Gaia DR3 reference epoch of J2016.0, and are propagated from the
# original values using the available proper motions. In the case of companions to higher-mass
# stars, sometimes the given coordinates are those of the primary - those have been offset by the
# separation and position angle of the companion, or RA and declination offsets, sourced from the
# notes column or the literature. Multiple systems whose components are listed separately in the
# main table and have distances that are divergent but consistent within 3-sigma use the mean of the
# individual measurements.

# Temperatures and radii, if not provided in the fundamental properties table, are interpolated from
# the BHAC15 (https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/) or COND03
# (https://perso.ens-lyon.fr/isabelle.baraffe/COND03_models) isochrones by Baraffe et al., using age
# and either luminosity or temperature, estimated from empirical relations. AbsMag is estimated the
# same way. Objects with a gravity classification of FLD-G have overestimated radii/underestimated
# temperatures given, so these have been recalculated.

# The relations of Sanghi et al. (2023), ApJ 959, 63 are used to estimate the luminosity of field
# dwarfs from absolute magnitude in different bands, if the parallax is precise to within 12.5%, or
# from the spectral type. Young objects (up to 300 Myr) use the corresponding absolute-magnitude
# relations, or the spectral type-temperature polynomial from that source. For cooler field dwarfs
# (Teff <~ 1000 K), the relation from W2-band magnitude by Leggett et al. (2025), ApJ 991, 193 is
# preferred. M dwarfs up to 10 Myr use the Teff scale from Herczeg and Hillenbrand (2014), ApJ 786,
# 97. Subdwarfs (up to type L7) use the Teff relation by Zhang et al. (2018), MNRAS 479, 1383; for
# objects of the sd subclass, an average of the field and subdwarf relations is taken. For the cases
# which fall outside the ranges of each of those relations, the field estimates are used, otherwise
# the Teff scale from the stellar properties table by Pecaut and Mamajek (2013)
# (http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt) is used.
"""

HEADER_BINS = """# Ultra-cool-dwarf binaries/multiples, based on the UltracoolSheet:
# https://zenodo.org/records/10573247

# Orbits are taken from the literature, unless there is no information on the masses or spectral
# types of the components. Published masses are used if they're directly measured (e.g. from orbital
# dynamics) or informed by some measured quantity, such as the system total mass or mass ratio.

# Systems without orbits take the coordinates from the main table and supplementary sources as those
# of the primary components, with the secondary positions being calculated from the separation and
# position angle, or from RA and declination offsets, sourced from the notes column of the main
# table or the literature.

# Magnitudes and radii are interpolated from the BHAC15 or COND03 isochrones the same way as for the
# main catalog. Missing masses are also estimated this way, and they were used for estimating the
# barycenter positions.
"""
