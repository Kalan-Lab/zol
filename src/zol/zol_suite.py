#!/usr/bin/env python3
"""
ZOL Suite Unified Interface

This module provides a unified command-line interface for all ZOL tools.
It allows users to run any of the main ZOL programs through a single entry point.
"""

import argparse
import sys
import subprocess
import os
from pathlib import Path
from typing import List
from rich_argparse import RawTextRichHelpFormatter, RawDescriptionRichHelpFormatter


# Custom descriptions for each program
PROGRAM_DESCRIPTIONS = {
    "zol": "Perform comparative gene cluster analysis.",
    "fai": "Find additional instances of gene clusters in a genome database using\n"
            "flexible alignment and synteny criteria.",
    "prepTG": "Prepare a database of target genomes for searches using fai.",
    "cgc": "Visualization of zol results along a consensus gene cluster sequence.",
    "cgcg": "Visualization of zol results as a graphical network.",
    "abon": "Automated analysis of conservation/novelty for a sample's biosynthetic\n"
            "gene clusters.",
    "apos": "Automated analysis of conservation/novelty for a sample's plasmids.",
    "atpoc": "Automated analysis of conservation/novelty for a sample's prophages.",
    "salt": "Support assessment for lateral transfer of gene clusters (experimental).",
    "zol-scape": "Run zol analysis on BiG-SCAPE results."
}


def get_zol_exec_path() -> Path:
    zol_exec_path = os.environ.get("ZOL_EXEC_PATH")
    if not zol_exec_path:
        print("Error: The environment variable ZOL_EXEC_PATH is not set. Please set it to the directory containing ZOL executables.")
        sys.exit(1)
    exec_path = Path(zol_exec_path)
    if not exec_path.exists() or not exec_path.is_dir():
        print(f"Error: ZOL_EXEC_PATH '{zol_exec_path}' does not exist or is not a directory.")
        sys.exit(1)
    return exec_path


def get_available_programs() -> List[str]:
    """Get list of available ZOL programs from ZOL_EXEC_PATH directory, limited to those in PROGRAM_DESCRIPTIONS."""
    bin_dir = get_zol_exec_path()
    all_programs = [f.name for f in bin_dir.iterdir() if f.is_file() and not f.name.startswith('.')]
    # Only include programs that are keys in PROGRAM_DESCRIPTIONS
    programs = [p for p in all_programs if p in PROGRAM_DESCRIPTIONS]
    
    # Force include 'zol' if it's in PROGRAM_DESCRIPTIONS, regardless of file status
    if 'zol' in PROGRAM_DESCRIPTIONS and 'zol' not in programs:
        programs.append('zol')
    
    return sorted(programs)


def run_program(program_name: str, args: List[str]) -> int:
    """
    Run a ZOL program with the given arguments, using the path from ZOL_EXEC_PATH.
    Args:
        program_name: Name of the program to run
        args: List of arguments to pass to the program
    Returns:
        Exit code from the program
    """
    bin_dir = get_zol_exec_path()
    program_path = bin_dir / program_name
    if not program_path.exists():
        print(f"Error: Program '{program_name}' not found in {bin_dir}.")
        return 1

    # Skip executability check for 'zol' to handle potential corruption issues
    if program_name != 'zol' and not os.access(program_path, os.X_OK):
        print(f"Error: Program '{program_name}' is not executable.")
        return 1
    try:
        cmd = [str(program_path)] + args
        result = subprocess.run(cmd, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {program_name}: {e}")
        return 1


def main() -> int:
    programs = get_available_programs()
    parser = argparse.ArgumentParser(
        description="""
The zol suite - a comprehensive bioinformatics toolkit for gene cluster analysis.

/================================\\
|| ________      ________       ||
|||\\_____  \\    |\\   ____\\      ||
|| \\|___/  /|   \\ \\  \\___|_     ||
||     /  / /    \\ \\_____  \\    ||
||    /  /_/__    \\|____|\\  \\   ||
||   |\\________\\    ____\\_\\  \\  ||
||    \\|_______|   |\\_________\\ ||
||                 \\|_________| ||
\\================================/

This interface provides access to all ZOL tools through a single command-line interface.
Each tool has its own specific arguments and functionality.

Typical order of operations: 

1.) Run [bold magenta]prepTG[/bold magenta] to prepare a database of target genomes for searches using fai.
2.) Run [bold magenta]fai[/bold magenta] to find additional instances of a gene cluster of interest in the prepTG database.
3.) Run [bold magenta]zol[/bold magenta] to perform comparative gene cluster analysis on the results from fai.
4.) Run [bold magenta]cgc[/bold magenta] and [bold magenta]cgcg[/bold magenta] to visualize the results from zol.

For help with a specific program, use: zol-suite <program> --help
        """,
        formatter_class=RawTextRichHelpFormatter,
    )
    parser.add_argument('--list-programs', action='store_true', help='List all available programs and exit')
    parser.add_argument('--version', action='version', version='ZOL Suite 1.6.0')
    subparsers = parser.add_subparsers(dest='program', metavar='<program>', help='ZOL program to run')

    # Only add subparsers for programs in PROGRAM_DESCRIPTIONS
    for prog in programs:
        desc = PROGRAM_DESCRIPTIONS.get(prog, f"Run {prog} (see `{prog} --help` for options)")
        subparsers.add_parser(
            prog,
            help=desc,
            description=desc,
            formatter_class=RawDescriptionRichHelpFormatter
        )

    # If zol-suite is invoked with no arguments, show help and exit 0
    if len(sys.argv) == 1:
        parser.print_help()
        return 0

    # Parse up to the subcommand, then pass the rest as args
    if len(programs) == 0:
        parser.print_help()
        print("\nError: No programs found in the directory specified by ZOL_EXEC_PATH.")
        return 1
    if len(sys.argv) > 1 and sys.argv[1] in programs:
        prog = sys.argv[1]
        prog_args = sys.argv[2:]
        return run_program(prog, prog_args)
    else:
        args = parser.parse_args()
        if args.list_programs:
            print("Available ZOL programs:")
            for program in programs:
                print(f"  {program}")
            return 0
        # If no subcommand or invalid, show help and available programs
        parser.print_help()
        print("\nAvailable ZOL programs:")
        for program in programs:
            print(f"  {program}")
        return 1


if __name__ == '__main__':
    sys.exit(main()) 