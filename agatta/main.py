# -*- coding: utf-8 -*-

"""
Agatta: Three-item analysis Python package
Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

Agatta is a set of tools in the cladistic framework to perform
three-item analysis and associated operations in cladistics in python.

https://github.com/vrineau/agatta

Usage:
    agatta analysis         <file>...        [options]
    agatta tripdec          <file>...        [options]
    agatta convert          <file>...        [options]
    agatta fp               <file>...        [options]
    agatta consensus        <file>           [options]
    agatta describetree     <file>           [options]
    agatta hmatrix          <file> ...       [options]
    agatta support          <file> <file>... [options]
    agatta chartest         <file> <file>... [options]
    agatta nri              <file> <file>... [options]
    agatta standardisation  <file> <file>... [options]
    agatta help             <command>
    agatta -h | --help
    agatta --version

Options:
    -h --help             Show general help information
    -s                    Silent mode
    -v                    Verbose mode
    --analysis=<type>     Type of searching method [default: heuristic]
    --chardec             Decompose trees into components
    --chartest            Compute hierarchical character states test
    --consensus=<type>    Specify type of consensus
    --directory=<dir>     Directory for output chartest files [default: ./]
    --detailed_tripdec    Compute a detailed table for triplet decomposition
    --filetype=<type>     Choose a tree/triplet/hmatrix file [default: trees]
    --index=<type>        Specify type of index to use [default: ri]
    --log                 Add a string to save an analysis log file
    --repetitions=<type>  Specify how to remove repetitions [default: TMS]
    --replicates=<int>    Number of replicates for an analysis [default: 1000]
    --parallel=<type>     Number of parallel cores to use [default: auto]
    --pdf                 Save a pdf file for chartest or nri results
    --prefix=<name>       Prefix of the output file [default: agatta_out]
    --ri                  Calculate Retention Index (Kitching et al. 1998)
    --rosetta=<path>      Path of the file for taxa conversion
    --showtaxanames       Print the list of terminals
    --softpath=<path>     Path of the software used
    --software=<type>     Software used for the pipeline [default: paup]
    --taxarep1=<path>     Table for leaves replacement for file 1
    --taxarep2=<path>     Table for leaves replacement for file 2
    --weighting=<type>    Specify the type of triplet weighting [default: FW]
    --rnri_codes=<path>   Compute NRI with a split of characters given in csv       
    --rnri_totaltree      Compute a pdf with absolute triplet support by node
    --rnri_rescaling      Rescale values in rnri_totaltree

"""

from .ini import helper
from .ini import hmatrix_several
from .ini import checkargs
from .ini import standardisation
from .ini import character_extraction
from .ini import wrapper_character_extraction
from .analysis import main_tripdec
from .analysis import del_replications_forest
from .interpret import RI_path
from .interpret import triplet_distance
from .interpret import rcc
from .interpret import constrict
from .interpret import describe_forest
from .interpret import chartest
from .interpret import NRI
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
            agatta_analysis(arguments["<file>"],  # one or several input files
                            arguments["--softpath"],
                            arguments["--software"],
                            arguments.get("--taxarep1", False),
                            arguments["--repetitions"],
                            arguments["--weighting"],
                            arguments["--parallel"],
                            arguments["--prefix"],
                            arguments["--analysis"],
                            arguments["--replicates"],
                            arguments.get("--rosetta", False),
                            arguments["--chartest"],
                            arguments["--ri"],
                            arguments.get("--consensus", False),
                            arguments["--pdf"],
                            arguments.get("--detailed_tripdec", False),
                            arguments.get("-v", False))

        # triplet decomposition
        elif arguments["tripdec"]:
            main_tripdec(arguments["<file>"],  # one or several input files
                         arguments["--prefix"],
                         arguments.get("--taxarep1", False),
                         arguments["--weighting"],
                         arguments["--parallel"],
                         arguments.get("--detailed_tripdec", False),
                         arguments["--repetitions"],
                         arguments.get("-v", False))

        # retention index calculation
        elif arguments["support"]:
            # retention index
            if arguments["--index"] == "ri":

                RI_path(arguments["<file>"][0],
                   arguments["<file>"][1:],  # one or several input files
                   arguments.get("--taxarep1", False),
                   arguments.get("--taxarep2", False),
                   arguments["--repetitions"],
                   arguments["--weighting"],
                   arguments["--prefix"]+".txt")

            # inter-tree retention index
            elif arguments["--index"] == "tripdistance":  # only two trees

                triplet_distance(list(character_extraction(
                    arguments["<file>"][0],
                    arguments.get("--taxarep1", False),
                    verbose=False).keys())[0],
                    list(character_extraction(
                        arguments["<file>"][1],
                        arguments.get("--taxarep2", False),
                        verbose=False).keys())[0],
                    arguments["--prefix"],
                    arguments["--weighting"])

        # character states testing procedure
        elif arguments["chartest"]:
            chartest(arguments["<file>"][0],
                     arguments["<file>"][1:],  # one or several input files
                     arguments.get("--taxarep1", False),
                     arguments.get("--taxarep2", False),
                     arguments["--repetitions"],
                     arguments["--prefix"],
                     arguments["--pdf"],
                     arguments.get("-v", False))

        # nodal retention index
        elif arguments["nri"]:
            NRI(arguments["<file>"][0],
                     arguments["<file>"][1:],  # one or several input files
                     arguments.get("--taxarep1", False),
                     arguments.get("--taxarep2", False),
                     arguments["--prefix"],
                     arguments.get("--rnri_codes", False),
                     arguments["--weighting"],
                     arguments["--repetitions"],
                     arguments.get("--rnri_totaltree", True),
                     arguments.get("--rescaling", False),
                     arguments["--pdf"])

        # convert
        elif arguments["convert"]:
            convert(arguments["<file>"],   # one or several input files
                    arguments["--filetype"],
                    arguments["--prefix"],
                    arguments["--parallel"],
                    arguments["--weighting"],
                    arguments["--analysis"],
                    arguments.get("--taxarep1", False),
                    arguments["--replicates"],
                    arguments.get("--log", False),
                    arguments["--software"],
                    False,
                    arguments["--repetitions"],
                    arguments.get("-v", False))

        # free-subtree paralogy analysis
        elif arguments["fp"]:
            del_replications_forest(wrapper_character_extraction(
                             arguments["<file>"],  # one or several input files
                             arguments.get("--taxarep1", False),
                             prefix=arguments["--prefix"],
                             verbose=arguments.get("-v", False)),
                             method=arguments["--repetitions"],
                             prefix=arguments["--prefix"],
                             verbose=arguments.get("-v", False))

        # consensus
        elif arguments["consensus"]:
            # reduced cladistic consensus
            if arguments["--consensus"] == "rcc":

                rcc(list(character_extraction(
                    arguments["<file>"][0],  # only one file
                    arguments.get("--taxarep1", False)).keys()),
                    arguments["--prefix"])

            # strict consensus
            else:
                constrict(list(character_extraction(
                          arguments["<file>"][0],   # only one file
                          arguments.get("--taxarep1", False)).keys()),
                          arguments["--prefix"])

        # describe trees
        elif arguments["describetree"]:   # only one file
            describe_forest(character_extraction(arguments["<file>"][0], 
                                                 info_tree=False),
                            arguments["--prefix"],
                            arguments.get("--showtaxanames", False))

        # standardisation
        elif arguments["standardisation"]:
            standardisation(arguments["<file>"][1:],  # one or several input
                              arguments["<file>"][0],
                              arguments["--prefix"],
                              verbose=arguments.get("-v", False))

        # transform hierarchical matrix into a tree list
        elif arguments["hmatrix"]:
            hmatrix_several(arguments["<file>"],  # one or several input
                    arguments["--prefix"],
                    arguments["--chardec"],
                    arguments["-v"])

        elif arguments["help"]:
            helper(arguments["<command>"])

        # display elapsed time
        elapsed_time = time.time() - start_time
        time_cptr = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print("elapsed time (total): {}".format(time_cptr))

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
