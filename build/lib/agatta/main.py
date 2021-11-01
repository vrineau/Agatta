# -*- coding: utf-8 -*-

"""
agatta: Three-item analysis Python package

agatta is a set of tools to to manipulate hierarchical characters and to 
perform three-item analysis and associated operations in python.

Usage:
    agatta analysis        <file>        [options]
    agatta tripdec         <file>        [options]
    agatta convert         <file>        [options]
    agatta fp              <file>        [options]
    agatta consensus       <file>        [options]
    agatta describetree    <file>        [options]
    agatta hmatrix         <file>        [options]
    agatta ri              <file> <file> [options]
    agatta chartest        <file> <file> [options]
    agatta standardisation <file> <file> [options]
    agatta help            <command>
    agatta -h | --help
    agatta --version

Options:
    -h --help            Show general help information
    -v                   Verbose mode
    -s                   Silent mode
    -i                   Replace taxa names for the 1st <file>
    -j                   Replace taxa names for the 2nd <file>
    --ri                 Calculate Retention Index (Kitching et al. 1998)
    --parallel=<type>    Number of parallel cores to use [default: auto]
    --prefix=<file>      Prefix of the output file [default: agatta_out]
    --directory=<dir>    Directory for output chartest files [default: ./]
    --software=<str>     Software used for the pipeline [default: tnt]
    --softpath=<path>    Path of the software used
    --weighting=<str>    Specify the type of triplet weighting [default: FW]
    --method=<str>       Specify the method to handle repetitions [default: VR]
    --consensus=<str>    Specify type of consensus [default: strict]
    --index=<str>        Specify type of index to use [default: ri]
    --chardec            Decompose trees into components
    --chartest           Compute hierarchical character states test (Cao 2007)
    --multiplier=<int>   Specifies a weight multiplicator [default: 1000000]
    --log                Add a string to save an analysis log file
    --filetype=<str>     Choose a tree file or a triplet file [default: trees]
    --showtaxanames      Print the list of terminals
    --analysis=<str>     Type of searching method [default: auto]
    --nrep=<int>         Number of replicates for an analysis [default: 1000]
    --rosette=<path>     Path of the file for taxa conversion
    --taxarep=<file>     Replacement file

"""

import sys
import os
import time
import docopt
from .ini import character_extraction
from .ini import phylo_to_areagram
from .ini import matrix_to_trees
from .ini import helper
from .analysis import del_replications_forest
from .analysis import main_tripdec
from .interpret import RI
from .interpret import ITRI
from .interpret import triplet_distance
from .interpret import character_states_test
from .interpret import constrict
from .interpret import rcc
from .interpret import describe_forest
from .out import agatta_analysis
from .out import convert
from .__version__ import __version__


def main():
    # docopt flags parsing
    arguments = docopt.docopt(__doc__, version=__version__)
    print(arguments)
    def coremain():

        start_time = time.time()

        print("agatta {}".center(80).format(__version__))
        print("Three-item analysis Python package".center(80))

        # Complete analysis
        if arguments["analysis"]:
            agatta_analysis(arguments["<file>"][0],
                            arguments["--softpath"],
                            arguments["--software"],
                            arguments.get("-i", False),
                            arguments["--method"],
                            arguments["--weighting"],
                            arguments["--parallel"],
                            arguments["--prefix"],
                            arguments["--analysis"],
                            arguments["--nrep"],
                            arguments.get("--rosette", False),
                            arguments["--chartest"],
                            arguments["--ri"],
                            arguments["--consensus"],
                            arguments.get("--pdf", False),
                            arguments.get("-v", False))

        # triplet decomposition
        elif arguments["tripdec"]:
            main_tripdec(arguments["<file>"][0],
                         arguments["--prefix"],
                         arguments.get("-i", False),
                         arguments["--weighting"],
                         arguments["--parallel"],
                         arguments.get("-v", False))

        # retention index calculation
        elif arguments["ri"]:
            # retention index
            if arguments["--index"] == "ri":

                RI(character_extraction(arguments["<file>"][0],
                                        arguments.get("-i", False)),
                   character_extraction(arguments["<file>"][1],
                                        arguments.get("-j", False)),
                   arguments["--weighting"],
                   output=arguments["--prefix"]+".txt")

            # inter-tree retention index
            elif arguments["--index"] == "itri":

                ITRI(list(character_extraction(
                    arguments["<file>"][0],
                    arguments.get("-i", False)).keys())[0],
                    list(character_extraction(
                        arguments["<file>"][1],
                        arguments.get("-j", False)).keys())[0],
                    arguments["--prefix"],
                    arguments["--weighting"])

            # triplet distance
            elif (arguments["--index"] == "itrisym_sum" or
                  arguments["--index"] == "itrisym_product"):

                triplet_distance(list(character_extraction(
                                arguments["<file>"][0],
                                arguments.get("-i", False)).keys())[0],
                                list(character_extraction(
                                arguments["<file>"][1],
                                arguments.get("-j", False)).keys())[0],
                                arguments["--prefix"],
                                arguments["--index"],
                                arguments["--weighting"])

        # character states testing procedure
        elif arguments["chartest"]:
            character_states_test(character_extraction(
                                  arguments["<file>"][0],
                                  arguments.get("-i", False)),
                                  character_extraction(
                                      arguments["<file>"][1],
                                      arguments.get("-j", False)),
                                  arguments["--prefix"],
                                  arguments.get("--pdf", False))

        # convert
        elif arguments["convert"]:
            convert(arguments["<file>"][0],
                    arguments["--filetype"],
                    arguments["--prefix"],
                    arguments["--parallel"],
                    arguments["--weighting"],
                    arguments["--analysis"],
                    arguments.get("-i", False),
                    arguments["--taxarep"],
                    arguments["--nrep"],
                    arguments.get("--log", False),
                    arguments["--software"],
                    arguments.get("-v", False),
                    arguments["--multiplier"])

        # free-subtree paralogy analysis
        elif arguments["fp"]:
            del_replications_forest(character_extraction(
                             arguments["<file>"][0],
                             arguments.get("-i", False)),
                             method=arguments["--method"],
                             prefix=arguments["--prefix"],
                             verbose=arguments["-v"])

        # consensus
        elif arguments["consensus"]:
            # strict consensus
            if arguments["--consensus"] == "strict":
                constrict(list(character_extraction(
                          arguments["<file>"][0],
                          arguments.get("-i", False)).keys()),
                          arguments["--prefix"])

            # reduced cladistic consensus
            elif arguments["--consensus"] == "rcc":

                rcc(list(character_extraction(
                    arguments["<file>"][0],
                    arguments.get("-i", False)).keys()),
                    arguments["--prefix"])

        # describe trees
        elif arguments["describetree"]:
            describe_forest(character_extraction(arguments["<file>"][0]),
                            arguments["--prefix"],
                            arguments.get("--showtaxanames", False))

        # standardisation
        elif arguments["standardisation"]:
            phylo_to_areagram(arguments["<file>"][0],
                              arguments["<file>"][1],
                              arguments["--prefix"],
                              verbose=arguments.get("-v", False))

        # transform hierarchical matrix into a tree list
        elif arguments["hmatrix"]:
            matrix_to_trees(arguments["<file>"][0],
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
