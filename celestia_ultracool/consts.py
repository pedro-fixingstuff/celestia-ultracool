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

# SpT-Teff scale for ~1 Myr old stars in IC 348 and Taurus (M1-M9) from Table 8
# of Luhman et al. (2003), 2003ApJ...593.1093L
SPT_TEFF_YOUNG = ((1, 2, 3, 4, 5, 6, 7, 8, 9),
                  (3705, 3560, 3415, 3270, 3125, 2990, 2880, 2710, 2400))

# Mean solar and Jovian radii
SOLAR_RADIUS = 695700.0
JUPITER_RADIUS = 69911.0

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

# Missing temperatures are estimated using the methods of Kirkpatrick et al. (2021), ApJS 253, 7 for
# field dwarfs. For young objects (up to 300 Myr), the relation from H-band absolute magnitude (from
# 2MASS or MKO systems) of Filippazzo et al. (2015), ApJ 810, 158 is used, and if there's no
# parallax precise to within 12.5% to calculate an absolute magnitude, the polynomial from Sanghi et
# al. (2023), ApJ 959, 63 is used. M dwarfs younger than 10 Myr (inclusive, for types earlier than
# M6) use the scale from Luhman et al. (2003), ApJ 593, 1093. Subdwarfs (up to type L7) use the
# relation by Zhang et al. (2018), MNRAS 479, 1383; for objects of the sd subclass, an average of
# the field and subdwarf relations is taken. For the cases which fall outside the ranges of those
# relations, the field estimates are used, otherwise the scale from the stellar properties table by
# Pecaut and Mamajek (2013) (http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt)
# is used.

# Radii (if not provided in the fundamental properties table) and magnitudes are interpolated from
# the BHAC15 (https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/) or COND03
# (https://perso.ens-lyon.fr/isabelle.baraffe/COND03_models) isochrones by Baraffe et al., using age
# and luminosity (if provided) or temperature. Objects with a gravity classification of FLD-G have
# overestimated radii/underestimated temperatures given, so these have been recalculated.
"""

HEADER_BINS = """# Ultra-cool-dwarf binaries/multiples, based on the UltracoolSheet:
# https://zenodo.org/records/10573247

# System coordinates from the main table and supplementary sources are taken as those of the primary
# components, with the secondary positions being calculated from the separation and position angle,
# or from RA and declination offsets, sourced from the notes column of the main table or the
# literature.

# Magnitudes and radii are interpolated from the BHAC15 or COND03 isochrones the same way as for the
# main catalog. Masses are also estimated this way, and they were used for estimating the barycenter
# positions.
"""
