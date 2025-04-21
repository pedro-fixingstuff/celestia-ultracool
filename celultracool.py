"""Builds ultra-cool dwarf catalogs for Celestia."""

import argparse
from celestia_ultracool.main import build_catalogs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Builds catalogs for Celestia from the UltracoolSheet.')
    parser.add_argument('-v', '--verbose', action='store_true', help="verbose output")
    parser.add_argument('-c', '--catalogs', action='store_true', help="write additional catalog designations (Pan-STARRS1, 2MASS, MKO, WISE)")
    parser.add_argument('-m', '--multiples', action='store_true', help="write companions and binary/multiple systems")
    parser.add_argument('-n', '--nearby', action='store_true', help="write objects within 10 pc but not included in nearstars.stc")
    args = parser.parse_args()
    build_catalogs(args.verbose, args.catalogs, args.multiples, args.nearby)
