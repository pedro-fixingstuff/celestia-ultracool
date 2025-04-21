"""Processes the UltracoolSheet into Celestia catalog files."""

import numpy as np
import pandas as pd
from pathlib import Path
import re

from . import consts
from .classes import System, Dwarf
from .utils import *


def build_catalogs(verbose: bool, write_catalogs: bool, write_multiples: bool, write_10pc: bool) -> None:
    """Builds Celestia .stc files for single and binary/multiple ultra-cool dwarfs."""
    home = Path.cwd()
    data_dir = home / 'data'
    output_dir = home / 'output'

    input_main = pd.read_csv(data_dir / 'UltracoolSheet - Main.csv')
    input_properties = pd.read_csv(data_dir / 'UltracoolSheet - FundamentalProperties.csv')

    input_suppl = pd.read_csv(data_dir / 'Supplemental_data.csv')
    input_exclusions = pd.read_csv(data_dir / 'Exclusions.csv')

    output = open(output_dir / 'ultracool.stc', 'w', encoding='utf-8')
    output.write(consts.HEADER)

    if write_multiples:
        input_binaries = pd.read_csv(data_dir / 'UltracoolSheet - Binaries.csv')
        input_triples = pd.read_csv(data_dir / 'UltracoolSheet - Triples+.csv')

        input_bin_suppl = pd.read_csv(data_dir / 'Binaries_supplemental_data.csv')

        output_bins = open(output_dir / 'ultracool_bins.stc', 'w', encoding='utf-8')
        output_bins.write(consts.HEADER_BINS)

    processed_count = 0

    for row in input_main.itertuples():
        if verbose:
            print('\nProcessing {}...'.format(row.name))

        # Check if object is in the exclusion list
        exclusion = input_exclusions[input_exclusions.name == row.name]
        if len(exclusion) > 0:
            if verbose:
                print('Object in exclusion list, skipping.')
            continue

        is_unresolved_multiple = row.multiplesystem_unresolved_in_this_table != 'N'
        is_multiple = (is_unresolved_multiple
                       or row.multiplesystem_resolved_in_this_table != 'N'
                       or row.has_higher_mass_companion != 'N')

        if not write_multiples and is_multiple:
            if verbose:
                print('Object is a binary/multiple system or a companion, skipping.')
            continue

        main_name = parse_name(row.name, is_multiple)
        main_name_stripped = strip_name(main_name)

        simbad_name = None

        suppl_row = input_suppl[input_suppl.name == row.name]
        has_suppl_data = len(suppl_row) != 0

        # Build name list
        names = [main_name]

        # Try to get the SIMBAD name from the supplemental table even if it's a "remove" value, thus
        # allowing the default SIMBAD name to be removed
        if has_suppl_data and pd.notna(suppl_row['name_simbad'].item()):
            simbad_name = parse_name(suppl_row['name_simbad'].item(), is_multiple)
            if simbad_name != 'remove' and strip_name(simbad_name) != main_name_stripped:
                names.append(simbad_name)

        # Otherwise, build the name list from scratch using the list of SIMBAD identifiers,
        # including Bayer, Flamsteed, variable-star and Gliese/GJ designations
        elif pd.notna(row.name_simbad):
            simbad_name = parse_name(row.name_simbad, is_multiple)

            simbad_identifiers = row.identifiers_simbad.split('|')

            bayer_flamsteed_designations = []
            variable_designation = None
            gj_designations = []

            for identifier in [str(id_) for id_ in simbad_identifiers]:
                if identifier.startswith('* '):
                    bayer_flamsteed_designations.append(parse_name(identifier))
                if identifier.startswith('V*') and not variable_designation:
                    variable_designation = parse_name(identifier)
                if identifier.startswith('GJ'):
                    # Check if two GJ designations are from different editions of the catalog
                    # (therefore having different IDs), as opposed to being the same, and referring
                    # to different components of a binary/multiple system
                    if not gj_designations or strip_name(gj_designations[-1]) != strip_name(identifier):
                        gj_designations.append(identifier)

            gj_designations = [parse_name(gj) for gj in gj_designations]

            names = gj_designations
            if variable_designation is not None:
                names.insert(0, variable_designation)
            names = bayer_flamsteed_designations + names

            if main_name_stripped not in [strip_name(name) for name in names]:
                names.append(main_name)

            if (strip_name(simbad_name).casefold() not in [strip_name(name).casefold() for name in names]):
                names.append(simbad_name)

        if has_suppl_data and pd.notna(suppl_row['aliases'].item()):
            aliases = suppl_row['aliases'].item()
            if aliases:
                names += aliases.split(':')

        # Store designations from large catalogs separately, to be appended in the case of single
        # objects or barycenters, while being left out for binary components
        catalogs = []

        if write_catalogs:
            if (pd.notna(row.designation_P1_formula) and str(row.designation_P1_formula).strip() not in names):
                catalogs.append(str(row.designation_P1_formula).strip())
            
            if (pd.notna(row.designation_2mass) and str(row.designation_2mass).strip() not in names):
                catalogs.append(str(row.designation_2mass).strip())
            
            if (pd.notna(row.designation_MKO) and str(row.designation_MKO).strip() not in names):
                catalogs.append(str(row.designation_MKO).strip())
            
            if (pd.notna(row.designation_WISE) and str(row.designation_WISE).strip() not in names):
                catalogs.append(str(row.designation_WISE).strip())

        if pd.notna(row.gucds_shortname):
            catalogs.append(row.gucds_shortname)

        offset_coords = False
        # Pick non-rounded coordinates from dataset-specific columns
        if row.source_j2000_formula in ['Gaia21a', 'Gaia23']:
            ra_epoch = row.ra_epoch_Gaia
            dec_epoch = row.dec_epoch_Gaia
            coord_decimal_digits = 7
            epoch = 2016.0

            # For companions that use their primary's coordinates, offset by their relative
            # astrometry (RA and dec. offsets, separation and position angle, or separation only)
            if row.astrom_Gaia == 'P':
                if has_suppl_data:
                    offset_coords = True
                    if pd.notna(suppl_row['offset_ew'].item()):
                        ra_epoch += suppl_row['offset_ew'].item() / np.cos(np.radians(dec_epoch)) / 3600
                    if pd.notna(suppl_row['offset_ns'].item()):
                        dec_epoch += suppl_row['offset_ns'].item() / 3600

                    elif pd.notna(suppl_row['sep'].item()):
                        sep = suppl_row['sep'].item()
                        pa = suppl_row['pa'].item()

                        ra_epoch += sep * np.sin(np.radians(pa)) / np.cos(np.radians(dec_epoch)) / 3600
                        dec_epoch += sep * np.cos(np.radians(pa)) / 3600

                elif pd.notna(row.sep_companion) and row.sep_companion:
                    offset_coords = True
                    dec_epoch += row.sep_companion / 3600
        elif row.source_j2000_formula == 'PS1':
            ra_epoch = row.ra_epoch_P1
            dec_epoch = row.dec_epoch_P1
            coord_decimal_digits = 4
            epoch = 2000 + (row.epoch_mjd_P1-51544.5) / 365.25
        elif row.source_j2000_formula == 'CatWISE':
            ra_epoch = row.ra_epoch_WISE
            dec_epoch = row.dec_epoch_WISE
            coord_decimal_digits = 7
            epoch = 2015.405
        elif row.source_j2000_formula == 'SIMBAD':
            ra_epoch = row.ra_j2000_simbad
            dec_epoch = row.dec_j2000_simbad
            coord_decimal_digits = 7
            epoch = 2000.0
        else:
            ra_epoch = row.ra_j2000_formula
            dec_epoch = row.dec_j2000_formula
            coord_decimal_digits = 4
            epoch = 2000.0

        # Pick proper motions supplemental table or from non-rounded dataset-specific columns
        if has_suppl_data and pd.notna(suppl_row['pmra'].item()):
            pm_ra = suppl_row['pmra'].item()
            pm_dec = suppl_row['pmdec'].item()
        if row.ref_pm_formula in ['Gaia21a', 'Gaia23']:
            pm_ra = row.pmra_Gaia
            pm_dec = row.pmdec_Gaia
        elif row.ref_pm_formula == 'Best20':
            pm_ra = row.pmra_UKIRT
            pm_dec = row.pmdec_UKIRT
        elif row.ref_pm_formula == 'Best18':
            pm_ra = row.pmra_P1
            pm_dec = row.pmdec_P1
        elif row.ref_pm_formula == 'Magn20':
            pm_ra = row.pmra_P1_PV34
            pm_dec = row.pmdec_P1_PV34
        elif row.ref_pm_formula == 'Maro21':
            pm_ra = row.pmra_catwise
            pm_dec = row.pmdec_catwise
        elif pd.notna(row.pmra_lit):
            pm_ra = row.pmra_lit
            pm_dec = row.pmdec_lit
        else:
            pm_ra = row.pmra_formula
            pm_dec = row.pmdec_formula

        ra = ra_epoch
        dec = dec_epoch
        # Propagate coordinates to the default epoch
        if pd.notna(pm_ra):
            ra += pm_ra / np.cos(np.radians(dec_epoch)) / 3600000 * (consts.DEFAULT_EPOCH - epoch)
        if pd.notna(pm_dec):
            dec += pm_dec / 3600000 * (consts.DEFAULT_EPOCH - epoch)

        use_abs_mag = False
        # Get parallax or distance from the supplemental or the main tables
        if has_suppl_data and pd.notna(suppl_row['plx'].item()):
            plx = suppl_row['plx'].item()
            plx_error = suppl_row['plxerr'].item()

            use_abs_mag = plx_error / plx <= 0.125

            dist_pc = 1000 / plx
            dist_note = 'from parallax ({})'.format(suppl_row['ref_plx'].item())
        elif has_suppl_data and pd.notna(suppl_row['dist'].item()):
            dist_pc = suppl_row['dist'].item()
            dist_note = suppl_row['dist_notes'].item()
        else:
            dist_pc, dist_error, dist_note = get_distance(row)

            # For components of resolved multiple systems with distance estimates that are divergent
            # but within 3-sigma from each other, use the mean of the individual measurements
            if row.multiplesystem_resolved_in_this_table != 'N':
                companion_name = row.multiplesystem_resolved_in_this_table.split(':')[0]

                companion_row = next(input_main[input_main.name == companion_name].itertuples())

                companion_dist, companion_dist_error, _ = get_distance(companion_row)

                if (companion_dist is not None and dist_pc != companion_dist
                    and abs(dist_pc - companion_dist) <= 3*np.sqrt(dist_error**2 + companion_dist_error**2)):
                    dist_pc = ((dist_pc / dist_error**2 + companion_dist / companion_dist_error**2)
                               / (1 / dist_error**2 + 1 / companion_dist_error**2))
                    dist_error = 1 / np.sqrt(1 / dist_error**2 + 1 / companion_dist_error**2)
                    dist_note += '; mean of system components'

            # Following Kirkpatrick et al. (2021), 2021ApJS..253....7K, if the parallax is known to
            # 12.5% or better, use the H-band absolute magnitude to calculate Teff
            use_abs_mag = pd.notna(row.plx_formula) and dist_error / dist_pc <= 0.125
        if dist_pc is not None:
            if not write_10pc and dist_pc <= 10:
                if verbose:
                    print('Object is within 10 pc, skipping.')
                continue
            dist = dist_pc * consts.PC_TO_LY
        else:
            dist = consts.PLACEHOLDER_DIST

        # Get spectral type (SpT) from the supplemental or the main table
        if has_suppl_data and pd.notna(suppl_row['sptnum'].item()):
            spt_num = suppl_row['sptnum'].item()
            if pd.notna(suppl_row['ref_spt'].item()):
                spt_note = suppl_row['ref_spt'].item()
            else:
                spt_note = None
        elif pd.notna(row.sptnum_formula):
            spt_num = row.sptnum_formula
            spt_note = None
        else:
            spt_num = consts.PLACEHOLDER_SPT_NUM
            spt_note = 'missing'

        is_subdwarf = spt_num < 0

        if has_suppl_data and pd.notna(suppl_row['age'].item()):
            age = suppl_row['age'].item()
            age_category = None
        else:
            age = row.age_singlevalue_gyr_formula
            age_category = row.age_category

        # Get metallicity subclass from the spectral type (optical or IR)
        if is_subdwarf:
            if ((pd.notna(row.spt_opt) and re.search(r'[eu]sd', row.spt_opt))
                or (pd.notna(row.spt_ir) and re.search(r'[eu]sd', row.spt_ir))):
                subclass = 'esd'
            else:
                subclass = 'sd'
        else:
            subclass = 'd'

        # Build InfoURL from SIMBAD name or use the SIMPLE link if provided
        if pd.notna(row.name_simbadable):
            infourl = ('https://simbad.cds.unistra.fr/simbad/sim-id?Ident='
                       + re.sub(r'\s+', '+', str(row.name_simbadable).replace('+', '%2B')))
        elif pd.notna(row.url_simpleDB):
            infourl = row.url_simpleDB
        else:
            infourl = None

        system = System(names, ra, dec, dist, spt_num, age, infourl, age_category, dist_note, spt_note)

        process_as_single = True
        # Process binary system components
        if is_unresolved_multiple:
            bin_suppl_row = input_bin_suppl[input_bin_suppl.name == row.name]
            has_bin_suppl_data = len(bin_suppl_row) != 0

            bin_query = input_binaries[input_binaries.name == row.name]
            if len(bin_query) > 0:
                for _ in range(1):  # stop if no separation is found
                    process_as_single = False

                    bin_row = next(bin_query.itertuples())

                    bin_ra_offset = 0.0
                    bin_dec_offset = 0.0
                    # Get coordinate offsets for the secondary
                    # Case 1: relative RA and dec. are given
                    if has_bin_suppl_data and pd.notna(bin_suppl_row['offset_ew_bin'].item()):
                        bin_ra_offset = bin_suppl_row['offset_ew_bin'].item() / np.cos(np.radians(dec)) / 3600000
                        bin_dec_offset = bin_suppl_row['offset_ns_bin'].item() / 3600000

                    # Case 2: both separation and position angle are given
                    elif has_bin_suppl_data and pd.notna(bin_suppl_row['pa_bin'].item()):
                        sep = bin_suppl_row['sep_bin'].item()
                        pa = bin_suppl_row['pa_bin'].item()

                        bin_ra_offset = sep * np.sin(np.radians(pa)) / np.cos(np.radians(dec)) / 3600000
                        bin_dec_offset = sep * np.cos(np.radians(pa)) / 3600000

                    # Case 3: only separation is given
                    elif has_bin_suppl_data and pd.notna(bin_suppl_row['sep_bin'].item()):
                        bin_dec_offset = bin_suppl_row['sep_bin'].item() / 3600000
                    elif pd.notna(bin_row.sep_bin):
                        bin_dec_offset = bin_row.sep_bin / 3600000

                    else:
                        output.write('\n# Binary/multiple, missing separation')
                        process_as_single = True
                        break

                    primary = Dwarf(system)
                    secondary = Dwarf(system)

                    # Use just "A" or "B" rather than "Aab" or "Bab" for barycenters
                    bin_designation = str(bin_row.designation_binary).replace('ab', '')

                    # Extract individual component designations
                    pri_designation = bin_row.designation_binary[:-1]
                    sec_designation = bin_row.designation_binary[:-2] + bin_row.designation_binary[-1]

                    bary_names = []
                    pri_names = []
                    sec_names = []
                    for name in names:
                        name = name.removesuffix(bin_row.designation_binary).removesuffix(bin_designation)
                        # Some primaries of wide (that is, resolved in the main table) binaries are
                        # themselves binary, so an 'A' letter means it's a companion to another
                        # object not addressed here and should not be removed
                        if pri_designation != 'A':
                            name = name.removesuffix(pri_designation)
                        name = name.rstrip()

                        # Only append the component letter to the barycenter designations if the
                        # binary is the secondary of a larger system
                        if not bin_designation.startswith('A'):
                            bary_names.append(name + ' ' + bin_designation)
                        else:
                            bary_names.append(name)

                        # Decompose multi-letter binary designations (e.g. "Bab")
                        if re.search(r'[aA-Z][bB-Z]$', name) is not None:
                            pri_name = name[:-1]
                            sec_name = name[:-2] + name[-1]
                        # If there's a component letter that's not been removed previously, append
                        # lowercase letters (e.g. "B" becomes "Ba"/"Bb")
                        elif re.search(r' [A-Z]$', name) is not None:
                            pri_name = name + 'a'
                            sec_name = name + 'b'
                        else:
                            pri_name = name + ' ' + pri_designation
                            sec_name = name + ' ' + sec_designation

                        pri_names.append(pri_name)
                        sec_names.append(sec_name)

                    system.names = bary_names + catalogs
                    primary.names = pri_names
                    secondary.names = sec_names

                    # Get component spectral types
                    if has_bin_suppl_data and pd.notna(bin_suppl_row['ref_spt'].item()):
                        primary.spt_num = bin_suppl_row['sptnum_pri'].item()
                        primary.spt_note = bin_suppl_row['ref_spt'].item()

                        secondary.spt_num = bin_suppl_row['sptnum_sec'].item()
                        secondary.spt_note = bin_suppl_row['ref_spt'].item()
                    else:
                        if pd.notna(bin_row.sptnum_pri):
                            primary.spt_num = bin_row.sptnum_pri
                            primary.spt_note = None
                        else:
                            primary.spt_note = 'missing, combined spectral type used'

                        if pd.notna(bin_row.sptnum_sec):
                            secondary.spt_num = bin_row.sptnum_sec
                            secondary.spt_note = None
                        else:
                            secondary.spt_note = 'missing, combined spectral type used'

                    # Calculate the H-band absolute magnitude for the components in order to try
                    # estimating their Teff
                    if use_abs_mag:
                        if pd.notna(bin_row.ref_H_MKO_bin_formula) and pd.notna(bin_row.ref_H_2MASS_bin_formula):
                            # Pick values with the lowest errors
                            if bin_row.Herr_MKO_pri_formula <= bin_row.Herr_2MASS_pri_formula:
                                pri_appmag_h = bin_row.H_MKO_pri_formula
                            else:
                                pri_appmag_h = bin_row.H_2MASS_pri_formula

                            if bin_row.Herr_MKO_sec_formula <= bin_row.Herr_2MASS_sec_formula:
                                sec_appmag_h = bin_row.H_MKO_sec_formula
                            else:
                                sec_appmag_h = bin_row.H_2MASS_sec_formula

                        elif pd.notna(bin_row.H_MKO_pri_formula):
                            pri_appmag_h = bin_row.H_MKO_pri_formula
                            sec_appmag_h = bin_row.H_MKO_sec_formula
                        else:
                            pri_appmag_h = bin_row.H_2MASS_pri_formula
                            sec_appmag_h = bin_row.H_2MASS_sec_formula

                        if pd.notna(pri_appmag_h):
                            pri_h_mag = app_to_abs_mag(pri_appmag_h, dist_pc)
                            primary.estimate_teff(pri_h_mag, subclass)
                        if pd.notna(sec_appmag_h):
                            sec_h_mag = app_to_abs_mag(sec_appmag_h, dist_pc)
                            secondary.estimate_teff(sec_h_mag, subclass)

                    if primary.teff is None:
                        primary.estimate_teff(subclass=subclass)
                    if secondary.teff is None:
                        secondary.estimate_teff(subclass=subclass)

                    primary.estimate_absmag()
                    primary.estimate_radius()
                    primary.estimate_mass()

                    secondary.estimate_absmag()
                    secondary.estimate_radius()
                    secondary.estimate_mass()

                    # Calculate secondary and barycenter positions from system coordinates (taken to
                    # be those of the primary) and relative astrometry, as well as mass ratio for
                    # the barycenter
                    secondary.ra += bin_ra_offset
                    secondary.dec += bin_dec_offset
                    system.ra += bin_ra_offset / (1 + primary.mass / secondary.mass)
                    system.dec += bin_dec_offset / (1 + primary.mass / secondary.mass)

                    system.components.append(primary)
                    system.components.append(secondary)

                    system.write(output_bins, is_subdwarf=is_subdwarf, offset_coords=offset_coords)

            # Process triple system components
            triple_query = input_triples[input_triples.name == row.name]
            if len(triple_query) > 0:
                process_as_single = False

                triple_row = next(triple_query.itertuples())

                primary = Dwarf(system)
                secondary = Dwarf(system)

                system.names += catalogs
                primary.names = [name + ' A' for name in names]
                secondary.names = [name + ' B' for name in names]

                # Get component spectral types
                if has_bin_suppl_data and pd.notna(bin_suppl_row['sptnum_pri'].item()):
                    primary.spt_num = bin_suppl_row['sptnum_pri'].item()
                    primary.spt_note = bin_suppl_row['ref_spt'].item()
                else:
                    primary.spt_num = triple_row.sptnum_1
                    primary.spt_note = None

                if has_bin_suppl_data and pd.notna(bin_suppl_row['sptnum_sec'].item()):
                    secondary.spt_num = bin_suppl_row['sptnum_sec'].item()
                    secondary.spt_note = bin_suppl_row['ref_spt'].item()
                else:
                    secondary.spt_num = triple_row.sptnum_2
                    secondary.spt_note = None

                primary.estimate_teff()
                secondary.estimate_teff()

                primary.estimate_absmag()
                primary.estimate_radius()
                primary.estimate_mass()

                secondary.estimate_absmag()
                secondary.estimate_radius()
                secondary.estimate_mass()

                if has_bin_suppl_data and pd.notna(bin_suppl_row['pa_bin'].item()):
                    sep = bin_suppl_row['sep_bin'].item()
                    pa = bin_suppl_row['pa_bin'].item()

                    bin_ra_offset = sep * np.sin(np.radians(pa)) / np.cos(np.radians(dec)) / 3600000
                    bin_dec_offset = sep * np.cos(np.radians(pa)) / 3600000
                else:
                    bin_ra_offset = 0.0
                    bin_dec_offset = triple_row.sep_21 / 3600000
                secondary.ra += bin_ra_offset
                secondary.dec += bin_dec_offset

                system.components.append(primary)

                if pd.notna(triple_row.sep_31):
                    subsystem = System(names, secondary.ra, secondary.dec, dist, secondary.spt_num, age, infourl, age_category)
                    tertiary = Dwarf(subsystem)

                    subsystem.names = [name + ' BC' for name in names]
                    tertiary.names = [name + ' C' for name in names]

                    if has_bin_suppl_data and pd.notna(bin_suppl_row['sptnum_ter'].item()):
                        tertiary.spt_num = bin_suppl_row['sptnum_ter'].item()
                        tertiary.spt_note = bin_suppl_row['ref_spt'].item()
                    else:
                        tertiary.spt_num = triple_row.sptnum_3
                        tertiary.spt_note = None

                    tertiary.estimate_teff()

                    tertiary.estimate_absmag()
                    tertiary.estimate_radius()
                    tertiary.estimate_mass()

                    if has_bin_suppl_data and pd.notna(bin_suppl_row['pa_tri'].item()):
                        sep = bin_suppl_row['sep_tri'].item()
                        pa = bin_suppl_row['pa_tri'].item()

                        triple_ra_offset = sep * np.sin(np.radians(pa)) / np.cos(np.radians(dec)) / 3600000
                        triple_dec_offset = sep * np.cos(np.radians(pa)) / 3600000
                    else:
                        triple_ra_offset = 0.0
                        triple_dec_offset = triple_row.sep_32 / 3600000
                    tertiary.ra += triple_ra_offset
                    tertiary.dec += triple_dec_offset

                    subsystem.ra += triple_ra_offset / (1 + secondary.mass / tertiary.mass)
                    subsystem.dec += triple_dec_offset / (1 + secondary.mass / tertiary.mass)

                    system.ra += (subsystem.ra - ra) / (1 + primary.mass / (secondary.mass + tertiary.mass))
                    system.dec += (subsystem.dec - dec) / (1 + primary.mass / (secondary.mass + tertiary.mass))

                    subsystem.components.append(secondary)
                    subsystem.components.append(tertiary)
                    system.components.append(subsystem)
                else:
                    system.ra += bin_ra_offset / (1 + primary.mass / secondary.mass)
                    system.dec += bin_dec_offset / (1 + primary.mass / secondary.mass)

                system.write(output_bins, is_subdwarf=is_subdwarf, offset_coords=offset_coords)

        if process_as_single:
            dwarf = Dwarf(system)
            dwarf.names = names + catalogs

            # Get physical parameters (radius and Teff) from table
            properties_query = input_properties[
                (input_properties.name == row.name) | (input_properties.name_simbadable == row.name_simbadable)]
            if len(properties_query) > 0:
                properties_row = next(properties_query.itertuples())

                lum = properties_row.log_lbol_lsun

                # Objects with a gravity classification of FLD-G, for some reason, have
                # overestimated radii and underestimated temperatures, so their parameters are
                # recalculated
                if properties_row.age_category in ['FLD-G', 'FLD-G?']:
                    dwarf.estimate_radius(lum)
                    dwarf.teff, dwarf.teff_note = calculate_teff(lum, dwarf.radius / consts.SOLAR_RADIUS)
                else:
                    dwarf.teff = properties_row.teff_evo
                    dwarf.radius = properties_row.radius_evo * consts.JUPITER_RADIUS

                dwarf.estimate_absmag(lum)

            # Otherwise, estimate Teff from empirical relations and use it to estimate parameters
            # from evolutionary models
            else:
                # Calculate the H-band absolute magnitude (from 2MASS or MKO photometric systems) to
                # estimate Teff. Unresolved systems are excluded, since they appear overluminous
                # compared to single objects of the same spectral type
                if use_abs_mag and not is_unresolved_multiple:
                    if pd.notna(row.Herr_MKO) and pd.notna(row.Herr_2MASS):
                        # Pick value with the lowest error
                        if row.Herr_MKO <= row.Herr_2MASS:
                            appmag_h = row.H_MKO
                        else:
                            appmag_h = row.H_2MASS
                    elif pd.notna(row.H_MKO):
                        appmag_h = row.H_MKO
                    else:
                        appmag_h = row.H_2MASS

                    if pd.notna(appmag_h):
                        h_mag = app_to_abs_mag(appmag_h, dist_pc)
                        dwarf.estimate_teff(h_mag, subclass)

                if dwarf.teff is None:
                    dwarf.estimate_teff(subclass=subclass)
                dwarf.estimate_absmag()
                dwarf.estimate_radius()

            dwarf.write(output, is_subdwarf=is_subdwarf, coord_decimal_digits=coord_decimal_digits, offset_coords=offset_coords)

        processed_count += 1
        if verbose:
            print('{} of {} object(s) processed.'.format(processed_count, len(input_main)))

    output.close()
    if write_multiples:
        output_bins.close()
