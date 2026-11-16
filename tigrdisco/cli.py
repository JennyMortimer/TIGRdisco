"""The tigrdisco command."""

import argparse

from . import __version__

COMMANDS = """commands:
  tigrdisco-fasta   find Tas genes in a nucleotide FASTA and search them for TIGR arrays
  tigrdisco-db      search GenBank records that carry a CENTRAL anchor
  tigrdisco-export  write Excel reports for a finished run
"""


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="tigrdisco",
        description="TIGRdisco: discovery of TIGR arrays and Tas proteins.",
        epilog=COMMANDS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"TIGRdisco {__version__}")
    parser.parse_args(argv)
    parser.print_help()


if __name__ == "__main__":
    main()
