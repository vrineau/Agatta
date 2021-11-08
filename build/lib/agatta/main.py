# -*- coding: utf-8 -*-

"""
Agatta: Three-item analysis Python package

Agatta is a set of tools in the cladistic framework to perform
three-item analysis and associated operations in cladistics in python.

https://github.com/vrineau/agatta

Usage:
    agatta analysis         <file>        [options]
    agatta tripdec          <file>        [options]
    agatta convert          <file>        [options]
    agatta fp               <file>        [options]
    agatta consensus        <file>        [options]
    agatta describetree     <file>        [options]
    agatta hmatrix          <file>        [options]
    agatta support          <file> <file> [options]
    agatta chartest         <file> <file> [options]
    agatta standardisation  <file> <file> [options]
    agatta help             <command>
    agatta -h | --help
    agatta --version

Options:
    -h --help            Show general help information
    -s                   Silent mode
    -v                   Verbose mode
    --analysis=<str>     Type of searching method [default: auto]
    --chardec            Decompose trees into components
    --chartest           Compute hierarchical character states test (Cao 2007)
    --consensus=<str>    Specify type of consensus [default: strict]
    --directory=<dir>    Directory for output chartest files [default: ./]
    --filetype=<str>     Choose a tree file or a triplet file [default: trees]
    --index=<str>        Specify type of index to use [default: ri]
    --log                Add a string to save an analysis log file
    --method=<str>       Specify how to remove repetitions [default: Rineau]
    --multiplier=<int>   Specifies a weight multiplicator [default: 1000000]
    --nrep=<int>         Number of replicates for an analysis [default: 1000]
    --parallel=<type>    Number of parallel cores to use [default: auto]
    --pdf=<path>         Specifies a path were to save pdf files for chartest
    --prefix=<file>      Prefix of the output file [default: agatta_out]
    --ri                 Calculate Retention Index (Kitching et al. 1998)
    --rosette=<path>     Path of the file for taxa conversion
    --showtaxanames      Print the list of terminals
    --softpath=<path>    Path of the software used
    --software=<str>     Software used for the pipeline [default: agatta]
    --taxarep1=<file>    Table for leaves replacement for file 1
    --taxarep2=<file>    Table for leaves replacement for file 2
    --weighting=<str>    Specify the type of triplet weighting [default: FW]

"""

from .ini import helper
from .ini import hmatrix
from .ini import checkargs
from .ini import standardisation
from .ini import character_extraction
from .analysis import main_tripdec
from .analysis import del_replications_forest
from .interpret import RI
from .interpret import ITRI
from .interpret import rcc
from .interpret import constrict
from .interpret import describe_forest
from .interpret import triplet_distance
from .interpret import chartest
from .out import agatta_analysis
from .out import convert
from .__version__ import __version__
import os
import sys
import time
import docopt


def main():
    # docopt flags parsing
    arguments = docopt.docopt(__doc__, version=__version__)

    def coremain():

        start_time = time.time()

        print("Agatta {}".center(80).format(__version__))
        print("Three-item analysis Python package".center(80))
        print()

        # check parsing
        checkargs(arguments)

        # Complete analysis
        if arguments["analysis"]:
            agatta_analysis(arguments["<file>"][0],
                            arguments["--softpath"],
                            arguments["--software"],
                            arguments.get("--taxarep1", False),
                            arguments["--method"],
                            arguments["--weighting"],
                            arguments["--parallel"],
                            arguments["--prefix"],
                            arguments["--analysis"],
                            arguments["--nrep"],
                            arguments.get("--rosette", False),
                            arguments["--chartest"],
                            arguments["--ri"],
                            arguments.get("--consensus", False),
                            arguments.get("--pdf", False),
                            arguments.get("-v", False))

        # triplet decomposition
        elif arguments["tripdec"]:
            main_tripdec(arguments["<file>"][0],
                         arguments["--prefix"],
                         arguments.get("--taxarep1", False),
                         arguments["--weighting"],
                         arguments["--parallel"],
                         arguments.get("-v", False))

        # retention index calculation
        elif arguments["support"]:
            # retention index
            if arguments["--index"] == "ri":

                RI(character_extraction(arguments["<file>"][0],
                                        arguments.get("--taxarep1", False),
                                        verbose=False),
                   character_extraction(arguments["<file>"][1],
                                        arguments.get("--taxarep2", False),
                                        verbose=False),
                   arguments["--weighting"],
                   arguments["--prefix"]+".txt")

            # inter-tree retention index
            elif arguments["--index"] == "itri":

                ITRI(list(character_extraction(
                    arguments["<file>"][0],
                    arguments.get("--taxarep1", False),
                    verbose=False).keys())[0],
                    list(character_extraction(
                        arguments["<file>"][1],
                        arguments.get("--taxarep2", False),
                        verbose=False).keys())[0],
                    arguments["--prefix"],
                    arguments["--weighting"])

            # triplet distance
            elif (arguments["--index"] == "itrisym_sum" or
                  arguments["--index"] == "itrisym_product"):

                triplet_distance(list(character_extraction(
                                arguments["<file>"][0],
                                arguments.get("--taxarep1", False),
                                verbose=False).keys())[0],
                                list(character_extraction(
                                arguments["<file>"][1],
                                arguments.get("--taxarep2", False),
                                verbose=False).keys())[0],
                                arguments["--prefix"],
                                arguments["--index"],
                                arguments["--weighting"])

        # character states testing procedure
        elif arguments["chartest"]:
            chartest(arguments["<file>"][0],
                     arguments["<file>"][1],
                     taxarep1=arguments.get("--taxarep1", False),
                     taxarep2=arguments.get("--taxarep2", False),
                     prefix=arguments["--prefix"],
                     pdf_files=arguments.get("--pdf", False))

        # convert
        elif arguments["convert"]:
            convert(arguments["<file>"][0],
                    arguments["--filetype"],
                    arguments["--prefix"],
                    arguments["--parallel"],
                    arguments["--weighting"],
                    arguments["--analysis"],
                    arguments.get("--taxarep1", False),
                    arguments["--nrep"],
                    arguments.get("--log", False),
                    arguments["--software"],
                    arguments.get("-v", False),
                    arguments["--multiplier"])

        # free-subtree paralogy analysis
        elif arguments["fp"]:
            del_replications_forest(character_extraction(
                             arguments["<file>"][0],
                             arguments.get("--taxarep1", False)),
                             method=arguments["--method"],
                             prefix=arguments["--prefix"],
                             verbose=arguments["-v"])

        # consensus
        elif arguments["consensus"]:
            # strict consensus
            if arguments.get("--consensus", "strict"):
                constrict(list(character_extraction(
                          arguments["<file>"][0],
                          arguments.get("--taxarep1", False)).keys()),
                          arguments["--prefix"])

            # reduced cladistic consensus
            elif arguments["--consensus"] == "rcc":

                rcc(list(character_extraction(
                    arguments["<file>"][0],
                    arguments.get("--taxarep1", False)).keys()),
                    arguments["--prefix"])

        # describe trees
        elif arguments["describetree"]:
            describe_forest(character_extraction(arguments["<file>"][0]),
                            arguments["--prefix"],
                            arguments.get("--showtaxanames", False))

        # standardisation
        elif arguments["standardisation"]:
            standardisation(arguments["<file>"][0],
                              arguments["<file>"][1],
                              arguments["--prefix"],
                              verbose=arguments.get("-v", False))

        # transform hierarchical matrix into a tree list
        elif arguments["hmatrix"]:
            hmatrix(arguments["<file>"][0],
                    arguments["--prefix"],
                    arguments["--chardec"],
                    arguments["-v"])

        elif arguments["help"]:
            exit(helper(arguments["<command>"]))

        # display elapsed time
        elapsed_time = time.time() - start_time
        time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    # class silent mode
    class HiddenPrints:
        def __enter__(self):
            self._original_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')

        def __exit__(self, exc_type, exc_val, exc_tb):
            sys.stdout.close()
            sys.stdout = self._original_stdout

    # choice silent mode or not
    if arguments["-s"]:
        with HiddenPrints():
            coremain()
    else:
        coremain()
