# Ultra-cool dwarf (a.k.a. brown dwarf) catalogs for Celestia

Celestia forum thread: https://celestiaproject.space/forum/viewtopic.php?f=23&t=23911

This project consists of Python scripts for generating catalog files from the [UltracoolSheet](http://bit.ly/UltracoolSheet) (UCS) for [Celestia](https://github.com/CelestiaProject/Celestia), a 3D space simulator.

Pre-built files are provided in the *Releases* subpage.

## Features

- Support for single ultra-cool dwarfs and binary/multiple systems
- Collects designations from several catalogs (including variable-star and GJ identifiers)
- Missing parameters filled by interpolation of the BHAC15 and COND03 isochrones
- Handles objects already in Celestia's default catalogs (via the `Exclusions` table)
- Allows for overriding/updating some parameters in the UCS (via the `Supplemental_data` and `Binaries_supplemental_data` tables)

## Setup

The scripts require Python 3.7 or higher to run.

1. Clone this repository or download it from GitHub.

2. Download the following files from the [Zenodo](https://zenodo.org/records/13993077) UCS repository, and save them to this repository's `data` directory:
    - `UltracoolSheet - Main.csv`
    - `UltracoolSheet - Binaries.csv`
    - `UltracoolSheet - Triples+.csv`
    - `UltracoolSheet - FundamentalProperties.csv`

3. Create a virtual environment (make sure the terminal is in this project's directory):
    ```bash
    python3 -m venv .venv
    ```

4. Activate the virtual environment:
    ```bash
    # Linux/macOS
    source .venv/bin/activate

    # Windows
    .venv\Scripts\activate
    ```

3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

   ```bash
   python celultracool.py [-v] [-c] [-m] [-n]
   ```

The following arguments are supported:

- `-v`, `--verbose`: Enable verbose output during generation.
- `-c`, `--catalogs`: Include additional catalog designations (Pan-STARRS1, 2MASS, MKO, WISE).
- `-m`, `--multiples`: Include companions, and binary/multiple systems.
- `-n`, `--nearby`: Include objects within 10 parsecs but not in Celestia's nearstars.stc catalog (which is limited to 25 ly).

### Output

These generated .stc files are saved in the `output` directory:

- `ultracool.stc`: Main catalog of ultra-cool dwarfs.
- `ultracool_bins.stc`: Catalog of binary and multiple systems (if `--multiples` is enabled)

## Acknowledgements

This work has benefitted from The UltracoolSheet at http://bit.ly/UltracoolSheet, maintained by Will Best, Trent Dupuy, Michael Liu, Aniket Sanghi, Rob Siverd, and Zhoujian Zhang, and developed from compilations by [Dupuy & Liu (2012)](http://adsabs.harvard.edu/abs/2012ApJS..201...19D), [Dupuy & Kraus (2013)](http://adsabs.harvard.edu/abs/2013Sci...341.1492D), [Deacon et al. (2014)](https://ui.adsabs.harvard.edu/abs/2014ApJ...792..119D), [Liu et al. (2016)](http://adsabs.harvard.edu/abs/2016ApJ...833...96L), [Best et al. (2018)](http://adsabs.harvard.edu/abs/2018ApJS..234....1B), [Best et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021AJ....161...42B), [Sanghi et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023ApJ...959...63S), and [Schneider et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023AJ....166..103S).	

This work has made use of the SIMBAD database, operated at CDS, Strasbourg, France ([Wenger et al. 2000](https://ui.adsabs.harvard.edu/abs/2000A%26AS..143....9W)).

## References

### Evolutionary model grids

- Baraffe et al. (2003), A&A 402, p. 701-712. ["Evolutionary models for cool brown dwarfs and extrasolar giant planets. The case of HD 209458"](https://ui.adsabs.harvard.edu/abs/2003A&A...402..701B) ([link to data](https://perso.ens-lyon.fr/isabelle.baraffe/COND03_models))
- Baraffe et al. (2015), A&A 577, A42. ["New evolutionary models for pre-main sequence and main sequence low-mass stars down to the hydrogen-burning limit"](https://ui.adsabs.harvard.edu/abs/2015A&A...577A..42B) ([link to data](https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/))

### Empirical relations

- Filippazzo et al. (2015), ApJ 810, 158. ["Fundamental Parameters and Spectral Energy Distributions of Young and Field Age Objects with Masses Spanning the Stellar to Planetary Regime"](https://ui.adsabs.harvard.edu/abs/2015ApJ...810..158F)
- Kirkpatrick et al. (2021), ApJS 253, 7. ["The Field Substellar Mass Function Based on the Full-sky 20 pc Census of 525 L, T, and Y Dwarfs"](https://ui.adsabs.harvard.edu/abs/2021ApJS..253....7K)
- Luhman et al. (2003), ApJ 593, p. 1093-1115. ["A Census of the Young Cluster IC 348"](https://ui.adsabs.harvard.edu/abs/2003ApJ...593.1093L)
- Pecaut and Mamajek (2013), ApJS 208, 9. ["Intrinsic Colors, Temperatures, and Bolometric Corrections of Pre-main-sequence Stars"](https://ui.adsabs.harvard.edu/abs/2013ApJS..208....9P) ([link to data](https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt))
- Sanghi et al. (2023), ApJ 959, 63. ["The Hawaii Infrared Parallax Program. VI. The Fundamental Properties of 1000+ Ultracool Dwarfs and Planetary-mass Objects Using Optical to Mid-infrared Spectral Energy Distributions and Comparison to BT-Settl and ATMO 2020 Model Atmospheres"](https://ui.adsabs.harvard.edu/abs/2023ApJ...959...63S)
- Zhang et al. (2018), MNRAS 479, p. 1383-1391. ["Primeval very low-mass stars and brown dwarfs - III. The halo transitional brown dwarfs"](https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.1383Z)

Values in the `Supplemental_data` and `Binaries_supplemental_data` tables come from a multitude of sources, cited by their bibcodes. Fields without a cited source have been extracted from the UCS main table.
