# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from .ini import hmatrix
from .ini import standardisation
from .ini import taxa_extraction
from .ini import taxa_to_numbers
from .ini import character_extraction
from .ini import wrapper_character_extraction
from .ini import taxa_triplet_extraction
from .analysis import main_tripdec
from .analysis import rep_detector
from .analysis import triplet_extraction
from .analysis import del_replications_forest
from .interpret import RI
from .interpret import rcc
from .interpret import constrict
from .interpret import character_states_test
from .search import search_pipeline
import os
import sys
import time
import datetime
import warnings
import treeswift
from ete3 import Tree


def triplet_nexus_file(triplet_dict, character_dict, weighting, analysis,
                       prefix, nrep=1000, logfile=False):
    """
    Outputs a three-item analysis file in nexus format readable by PAUP*4 from
    a dictionary of triplets and their corresponding weights.

    Swofford, D. L. 2003. PAUP*. Phylogenetic Analysis Using Parsimony
    (*and Other Methods). Version 4. Sinauer Associates, Sunderland,
    Massachusetts.

    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    weighting : str
        Weighting scheme to use between FW, FWNL, MW, AW, NW.
    analysis : str
        Add a line in the nexus file to state the type of analysis between
        branch and bound ("bandb") or heuristic ("heuristic").
    prefix : str
        Prefix of the saving nexus file prefix.nex. The complete path can be
        used.
    nrep : int, optional
        Add a line in the nexus file to state the number of replicates in case
        of heuristic search. The default is 1000.
    logfile : bool, optional
        IF true, Add a line in the nexus file to compute a log file during the
        PAUP4 analysis. The default is False.

    Returns
    -------
    None.

    """

    # weights
    count_zeroweights = 0  # zero-weight triplets (round)

    if weighting != "NW":
        i = 1  # triplet identifier
        triplet_list_FW = str()  # weight lines
        deltriplets = dict()

        for trip, FW in triplet_dict.items():

            # weights management using maximum weight in PAUP4
            if weighting == "FW" or weighting == "FWNL" or weighting == "MW":
                if character_dict:
                    FWdec = len(str(len(character_dict)))
                else:
                    FWdec = 2
                nexw = round((float(FW)), 8 - FWdec)
                if nexw == 0:
                    deltriplets[trip] = 0  # if zero-weight
                    count_zeroweights += 1
                elif nexw.is_integer():
                    triplet_list_FW += str(int(FW))
                else:
                    triplet_list_FW += str(nexw)
            else:
                nexw = 1
                triplet_list_FW += str(int(FW))  # NW, AW, UW (integers)

            if nexw != 0:
                triplet_list_FW += ":"
                triplet_list_FW += str(i)
                triplet_list_FW += ", "
                i += 1

        for trip in deltriplets.keys():
            triplet_dict.pop(trip)

    if character_dict:
        taxa_list = taxa_extraction(character_dict)
    else:
        taxa_list = taxa_triplet_extraction(triplet_dict)

    # if not prefix:
    #     root = Tk()
    #     prefix = filedialog.asksaveasfilename(
    #         title="Save three-taxon analysis PAUP nexus file", filetypes=((
    #             "nexus files", "*.nex"), ("all files", "*.*")))
    #     root.withdraw()

    # Nexus file computation
    with open(prefix+".nex", "w") as nexus_file:
        nexus_file.write("#NEXUS")
        nexus_file.write("\nbegin data;")
        nexus_file.write("\nDimensions ntax={} nchar={};".format(
            str(1 + len(taxa_list)), str(len(triplet_dict))))
        nexus_file.write("\nFormat symbols=\"0 1\" missing=?;")
        nexus_file.write("\n")
        nexus_file.write("\nMatrix")
        nexus_file.write("\n")

        # matrix computation
        for taxa_name in taxa_list:
            nexus_file.write("\n"+taxa_name+" ")
            newline = ""

            # for each triplet
            for triplet_column in triplet_dict.keys():
                if (taxa_name not in
                        triplet_column.out_taxa | triplet_column.in_taxa):
                    newline += "?"
                elif taxa_name in triplet_column.in_taxa:
                    newline += "1"
                else:
                    newline += "0"

            nexus_file.write(newline)

        nexus_file.write("\nroot "+"0"*len(triplet_dict))  # all 0 root

        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")
        nexus_file.write("\n")
        nexus_file.write("\nBegin Paup;")

        # weights lines
        if weighting != "NW":
            nexus_file.write("\nwts")
            nexus_file.write("\n"+triplet_list_FW[:-2])

        nexus_file.write("\n;")
        nexus_file.write("\noutgroup root /only;")

        if logfile:
            nexus_file.write("\nlog file={}_paup.log;".format(prefix))

        nexus_file.write("\nset maxtrees = 1000 increase=auto;")

        if analysis == "bandb":
            nexus_file.write("\nbandb;")
        elif analysis == "heuristic":
            nexus_file.write("\nhsearch addseq=random nreps={};".format(nrep))

        nexus_file.write("\nroottrees;")
        nexus_file.write("\nsavetrees /file={}.tre format=newick;".format(
            prefix))

        if logfile:
            nexus_file.write("\nlog stop;")  # end of log save

        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")

        if count_zeroweights > 0:
            print("WARNING, " + count_zeroweights +
                  " zero weight triplets deleted due to rounding")


def triplet_tnt_file(triplet_dict, character_dict, weighting, analysis,
                     prefix, nrep=1000, logfile=False):
    """
    Outputs a three-item analysis file in TNT format (Goloboff et al. 2008)
    from a dictionary of triplets and their corresponding weights.

    Goloboff, P. A., Farris, J. S., & Nixon, K. C. (2008). TNT, a free program
    for phylogenetic analysis. Cladistics, 24(5), 774-786.

    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    weighting : str
        Weighting scheme to use between FW, FWNL, MW, AW, NW.
    analysis : str
        Add a line in the TNT file to state the type of analysis between
        branch and bound ("bandb") or heuristic ("heuristic").
    prefix : str
        Prefix of the saving TNT file prefix.tnt. The complete path can be
        used.
    nrep : int, optional
        Add a line in the TNT file to state the number of replicates in case
        of heuristic search. The default is 1000.
    logfile : bool, optional
        IF true, Add a line in the TNT file to compute a log file during the
        TNT analysis. The default is False.

    Returns
    -------
    None.

    """

    count_zeroweights = 0  # zero-weight triplets (round)

    if weighting != "NW":
        triplet_list_FW = "ccode "
        wfactor = 1000 / float(max(triplet_dict.values()))
        itrip = 1
        deltriplets = dict()

        # matrix computation
        if (weighting == "FW" or weighting == "FWNL" or weighting == "MW" or
                max(triplet_dict.values()) > 1000):  # fractional weights
            for triplet_column, FW in triplet_dict.items():  # for each triplet

                nexw = int(float(FW)*wfactor)
                if nexw == 0:
                    deltriplets[triplet_column] = 0  # if zero-weight triplet
                    count_zeroweights += 1
                else:
                    triplet_list_FW += " /{} ".format(str(nexw))
                    triplet_list_FW += str(itrip)

                    itrip += 1

        # if rescale weight not needed
        else:
            for triplet_column, FW in triplet_dict.items():  # for each triplet

                nexw = 1
                triplet_list_FW += " /{} ".format(str(FW))
                triplet_list_FW += str(itrip)

                itrip += 1

        triplet_list_FW += "*;"  # end of the ccode line
        for trip in deltriplets.keys():
            triplet_dict.pop(trip)

    # computation of the TNT file
    # if prefix is None:
    #     root = Tk()
    #     prefix = filedialog.asksaveasfilename(
    #         title="Save three-taxon analysis PAUP nexus file",
    #         filetypes=(("nexus files", "*.nex"), ("all files", "*.*")))
    #     root.withdraw()

    if character_dict:
        taxa_list = taxa_extraction(character_dict)
    else:
        taxa_list = taxa_triplet_extraction(triplet_dict)

    tntstring = "xread\n{} {}".format(
        str(len(triplet_dict)), (str(1 + len(taxa_list))))
    tntstring += "\n"

    # matrix computation
    for taxa_name in taxa_list:
        tntstring += "\n"+taxa_name+" "
        newline = ""

        # for each triplet
        for triplet_column in triplet_dict.keys():

            if (taxa_name not in
                    triplet_column.out_taxa | triplet_column.in_taxa):
                newline += "?"
            elif taxa_name in triplet_column.in_taxa:
                newline += "1"
            else:
                newline += "0"

        tntstring += newline

    tntstring += "\nroot "+"0"*len(triplet_dict)  # all zeros root line

    tntstring += "\n;"
    tntstring += "\n"

    # weights line
    if weighting != "NW":
        tntstring += "\n"+triplet_list_FW

    tntstring += "\n;"
    tntstring += "\n"

    if logfile:
        tntstring += "\nlog {}.tnt.log;".format(prefix)  # if log

    tntstring += "\ntaxname =;"

    if analysis == "bandb":
        tntstring += "\nienum;"

    elif analysis == "heuristic":
        tntstring += "\nhold 10000;"  # max trees
        tntstring += "\nmult=replic {};".format(str(nrep))  # heuristic
        tntstring += "\nbbreak=tbr;"  # swapping method
        tntstring += "\ncollapse [; collapse 1;"  # collapse null branches

    tntstring += "\nexport - {}.tre;".format(prefix)  # export results to file
    tntstring += "\nquit"

    if count_zeroweights > 0:
        print("WARNING: " + str(count_zeroweights) +
              " zero weight triplets deleted due to rounding")

    with open(prefix+".tnt", "w") as tnt_file:
        tnt_file.write(tntstring)


def triplet_tmc_file(triplet_dict, character_dict, prefix,
                     weighting="FW"):
    """
    Outputs a three-item analysis file readable by Triplet MaxCut
    (Sevillya et al., 2016) from a dictionary of triplets and their
    corresponding weights.

    Sevillya, G., Frenkel, Z., & Snir, S. (2016). Triplet MaxCut: a new toolkit
    for rooted supertree. Methods in Ecology and Evolution, 7(11), 1359-1365.

    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    prefix : str
        Prefix of the saving TNT file prefix.tnt. The complete path can be
        used.
    weighting : str
        Weighting scheme to use between FW, FWNL, MW, AW, NW.

    Returns
    -------
    None.

    """

    # if prefix is None:
    #     root = Tk()
    #     prefix = filedialog.asksaveasfilename(
    #         title="Save three-taxon analysis TMC file", filetypes=((
    #             "text files", "*.txt"), ("all files", "*.*")))
    #     root.withdraw()

    # file for translating taxa names to numbers
    taxa_list = taxa_extraction(character_dict)
    taxa_dict, taxa_convers = taxa_to_numbers(taxa_list)

    with open(prefix+"_taxa.txt", "w") as TMC_file:
        for convers_line in taxa_convers:
            TMC_file.write(convers_line+"\n")
        TMC_file.write("\n")

    # string file
    tmcstring = ''
    for triplet1 in triplet_dict.keys():

        in_tax_list = list(triplet1.in_taxa)
        (out_taxa_str,) = triplet1.out_taxa

        # if weights are fractions
        if weighting == "FW" or weighting == "FWNL" or weighting == "MW":
            if character_dict:
                FWdec = len(str(len(character_dict)))
            else:
                FWdec = 2

            tmcstr = taxa_dict[in_tax_list[0]] + ","
            + taxa_dict[in_tax_list[1]] + "|" + taxa_dict[out_taxa_str]
            + "[:" + str(round((float(triplet1.FW)), 8 - FWdec)) + "] "
            tmcstring += tmcstr

        # if weights are integers
        else:
            tmcstring += taxa_dict[in_tax_list[0]] + ","
            + taxa_dict[in_tax_list[1]] + "|" + taxa_dict[out_taxa_str]
            + "[:" + str(int(triplet1.FW)) + "] "

    # write TMC file
    with open(prefix+".tmc", "w") as TMC_file:
        TMC_file.write(tmcstring)


def triplet_wqfm_wtqmc_file(triplet_dict, prefix, weighting="FW", precision=7):
    """
    Outputs a three-item analysis file readable by wqfm (Mahbub et al., 2021)
    and wTREE-QMC (Han & Molloy, 2024) from a dictionary of triplets and their 
    corresponding weights. To allow a three-item analysis with quartets, 
    each quartet possesses a leaf which is the root.

    Han and Molloy, 2024, Improved robustness to gene tree incompleteness, 
    estimation errors, and systematic homology errors with weighted TREE-QMC, 
    bioRxiv, https://doi.org/10.1101/2024.09.27.615467.

    Mahbub, M., Wahab, Z., Reaz, R., Rahman, M. S., & Bayzid, M. (2021).
    wQFM: Highly Accurate Genome-scale Species Tree Estimation from Weighted
    Quartets. Bioinformatics.
    
    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).
    prefix : str, optional
        Prefix of the saving wqfm file.
    weighting : str, optional
        Weighting scheme to use between FW, FWNL, MW, AW, NW.
        The default is "FW".
    precision : int, optional
        Number of decimals to keep when rounding fractional weights.
        Increase this value for better precision.

    Returns
    -------
    None.

    """

    # core wqfm or wTREE-QMC file
    with open(prefix, "w") as wqfm_file:
        for triplet1, FW in triplet_dict.items():

            # if weights are fractions
            if weighting == "FW" or weighting == "FWNL" or weighting == "MW":
                wqfm_file.write("((" + list(triplet1.out_taxa)[0]
                                + ",root),(" + list(triplet1.in_taxa)[0] + ","
                                + list(triplet1.in_taxa)[1] + ")); "
                                + str(round(float(FW),precision))  # precision
                                + "\n")

            # if weights are integers
            else:
                wqfm_file.write("((" + list(triplet1.out_taxa)[0]
                                + ",root),(" + list(triplet1.in_taxa)[0]
                                + "," + list(triplet1.in_taxa)[1] + ")); "
                                + str(round(FW)) + "\n")


def triplet_to_file(triplet_dict, character_dict, prefix, analysis="heuristic",
                    nrep=1000, logfile=True, software="paup", weighting="FW"):
    """
    Outputs a three-item analysis file readable by TNT (Goloboff et al. 2008),
    by PAUP* (Swofford, 2003), by WQFM (Mahbub et al. 2021), or by wTREE-QMC 
    (Han & Molloy, 2024) from a dictionary of triplets and their corresponding 
    weights.

    Goloboff, P. A., Farris, J. S., & Nixon, K. C. (2008). TNT, a free program
    for phylogenetic analysis. Cladistics, 24(5), 774-786.

    Han and Molloy, 2024, Improved robustness to gene tree incompleteness, 
    estimation errors, and systematic homology errors with weighted TREE-QMC, 
    bioRxiv, https://doi.org/10.1101/2024.09.27.615467.
    
    Mahbub, M., Wahab, Z., Reaz, R., Rahman, M. S., & Bayzid, M. (2021).
    wQFM: Highly Accurate Genome-scale Species Tree Estimation from Weighted
    Quartets. Bioinformatics.

    Swofford, D. L. 2003. PAUP*. Phylogenetic Analysis Using Parsimony
    (*and Other Methods). Version 4. Sinauer Associates, Sunderland,
    Massachusetts.

    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. If set to none, a window appears to choose the location and name
        of the prefix for the nexus file. The default is None.
    analysis : str
        Add a line in the nexus file to state the type of analysis between
        branch and bound ("bandb") or heuristic ("heuristic").
    nrep : int, optional
        Add a line in the file to state the number of replicates in case
        of heuristic search. The default is 1000.
    logfile : bool, optional
        IF true, Add a line to outputs a log file during the  analysis.
        The default is False.
    software : str, optional
        Choose the software to use between 'paup' and 'tnt'.
        The default is "paup".
    weighting : str
        Weighting scheme to use between FW, FWNL, MW, AW, NW.
    Returns
    -------
    None.

    """

    print("Computing {} file".format(software))

    start = time.time()

    if software == "paup":
        triplet_nexus_file(triplet_dict, character_dict, weighting, analysis,
                           prefix, nrep, logfile)

    elif software == "tnt":
        triplet_tnt_file(triplet_dict, character_dict, weighting, analysis,
                         prefix, nrep, logfile)

    # elif software == "tmc":
    #     triplet_tmc_file(triplet_dict, character_dict, prefix, weighting)

    elif software == "wqfm":
        triplet_wqfm_wtqmc_file(triplet_dict, prefix+".wqfm", weighting)
        
    elif software == "wtree-qmc":
        triplet_wqfm_wtqmc_file(triplet_dict, prefix+".wtqmc", weighting)

    end = time.time()
    time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print("elapsed time (output file computation): {}".format(time_cptr))


def convert(infilelist, infiletype, prefix, parallel="auto", weighting="FW",
            analysis="heuristic", taxa_replacement=False, nreplicates=1000,
            logfile=True, software="paup", dec_detail=False, method="TMS", 
            verbose=True):
    """
    Outputs a three-item analysis file readable by TNT (Goloboff et al. 2008),
    by PAUP* (Swofford, 2003), by WQFM (Mahbub et al. 2021), or by
    wTREE-QMC (Han & Molloy, 2024) from :

        * a dictionary of triplets and their corresponding
          weights, e.g., computed from main_tripdec.
        * a dictionary of rooted trees.
        * a path to a text file containing newick trees.

    Goloboff, P. A., Farris, J. S., & Nixon, K. C. (2008). TNT, a free program
    for phylogenetic analysis. Cladistics, 24(5), 774-786.
    
    Han and Molloy, 2024, Improved robustness to gene tree incompleteness, 
    estimation errors, and systematic homology errors with weighted TREE-QMC, 
    bioRxiv, https://doi.org/10.1101/2024.09.27.615467.

    Mahbub, M., Wahab, Z., Reaz, R., Rahman, M. S., & Bayzid, M. (2021).
    wQFM: Highly Accurate Genome-scale Species Tree Estimation from Weighted
    Quartets. Bioinformatics.

    Swofford, D. L. 2003. PAUP*. Phylogenetic Analysis Using Parsimony
    (*and Other Methods). Version 4. Sinauer Associates, Sunderland,
    Massachusetts.

    Parameters
    ----------
    infilelist : dict or list of str
        Dictionary of ete3. Tree in keys or dictionary of triplets in keys and
        weights in values or list of paths of file(s) containing newick trees.
    infiletype : str
        Type of infile argument. Can be 'triplets' if infile is a dictionary
        of triplets or 'trees' if infile is a path to a tree file or a
        dictionary of trees.
    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. If set to none, a window appears to choose the location and name
        of the prefix for the nexus file. The default is None.
    parallel : str or int
        Option for choosing if the analysis is made using multiprocessing or
        not. This argument can be:

            * "not" if the user does not wan to use multiprocessing.
            * "auto" for automatic detection of the number of cpu available.
            * any integer corresponding to the number of cpu the user wants to
              allow to the analysis.

        This argument is only used if infiletype is 'trees'.
    weighting : str
        Weighting scheme to use between FW, FWNL, MW, AW, NW.
    analysis : str
        Add a line in the nexus file to state the type of analysis between
        branch and bound ("bandb") or heuristic ("heuristic").
    taxa_replacement : str, optional
        Path of a table file containing two columns. The first column
        corresponds to the names of the terminals of the newick stored in
        infile, and the second column corresponds to their names the user wants
        to obtain at the end. All separators accepted. Example:
             AA Diceras
             AB Valletia
             AC Monopleura
        The default is False (no replacement).
    nrep : int, optional
        Add a line in the file to state the number of replicates in case
        of heuristic search. The default is 1000.
    logfile : bool, optional
        IF true, Add a line to outputs a log file during the  analysis.
        The default is False.
    software : str, optional
        Choose the software to use between 'paup', 'tnt', 'wqfm', 'wtree-qmc'.
        The default is "paup".
    dec_detail : bool, optional
        If true, save a detailed table in csv format named 
        prefix.triplet_table.csv of triplet weights per character in addition 
        to the .triplet file. The default is False.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    verbose : bool, optional
        Verbose mode if True. The default is True.
    Returns
    -------
    None.

    """

    if infiletype == "triplets":
        triplet_dict = triplet_extraction(infilelist, taxa_replacement)
        character_dict = False

    elif infiletype == "trees":
        if type(infilelist) == list:
            character_dict = wrapper_character_extraction(infilelist, 
                                                          taxa_replacement,
                                                          prefix,
                                                          verbose)

        else:
            character_dict = infilelist

        triplet_dict = main_tripdec(infilelist, prefix, taxa_replacement,
                                    weighting, parallel, dec_detail, 
                                    method, verbose)

    triplet_to_file(triplet_dict, character_dict, prefix, analysis,
                    nreplicates, logfile, software, weighting)


def agatta_analysis(infilelist, software_path, software="paup",
                    taxa_replacement=False, method="TMS", weighting="FW",
                    parallel="auto", prefix="agatta_out", analysis="heuristic",
                    nrep=1000, rosetta=False, chartest=False, ri=False,
                    consensus=False, pdf_file=False, dec_detail=False, 
                    verbose=True):
    """
    Main function of the Agatta python package. Allow to perform a three-item
    analysis (Nelson & Platnick, 1991), e.g., in the context of systematics
    phylogenetics, using hierarchical characters (Cao et al. 2007), or in
    cladistic biogeography.

    The analysis can be performed using a text file containing a hierarchical
    matrix (see https://vrineau.github.io/AgattaDocs/Input%20files.html for 
    informations about the format) or newick rooted trees.

    Repetitions (polymorphism) are automatically removed according to:

        * The original algorithm of Nelson and Ladiges (1996) designed in the
          paradigm of cladistic biogeography.

        * The algorithm of Rineau et al. (in prep) designed to construct
          subtrees without repeated leaves while minimising the loss of
          information in terms of triplets. This is the default mode.

    The following options are available for analysing the results:

        * Computation of a consensus tree (strict consensus or reduced
          cladistic consensus; Wilkinson 1994)

        * Character state test procedure for hierarchical characters (Cao 2008
          modified by Rineau (2017)
          to analyse if character states are synapomorphies or homoplasies.

        * Retention index computation for all characters
          (Kitching et al., 1998).

          Cao, N. (2008). Analyse à trois éléments et anatomie du bois des
          Fagales Engl (Doctoral dissertation, Paris, Muséum national
          d'histoire naturelle).

          Cao, N., Zaragüeta Bagils, R., & Vignes-Lebbe, R. (2007).
          Hierarchical representation of hypotheses of homology.
          Geodiversitas, 29(1), 5-15.

          Kitching, I. J., Forey, P., Humphries, C., & Williams, D. (1998).
          Cladistics: the theory and practice of parsimony analysis.
          Oxford University Press.

          Nelson, G. J., & Ladiges, P. Y. (1996). Paralogy in cladistic
          biogeography and analysis of paralogy-free subtrees.
          American Museum novitates 3167.

          Nelson, G., & Platnick, N. I. (1991). Three‐taxon statements:
          a more precise use of parsimony?. Cladistics, 7(4), 351-366.

          Rineau, V. (2017). Un nouveau regard cladistique sur l'anatomie
          comparée, la phylogénie, la systématique et la paléoécologie des
          rudistes (Bivalvia, Hippuritida) (Doctoral dissertation,
          Université Pierre et Marie Curie).

          Rineau, V., Moncel, M-H., & Zeitoun, V. (in prep). Revealing
          evolutionary patterns behind homogeneity: the case of the Paleolithic
          assemblages from Notarchirico (Southern Italy).

          Wilkinson, M. (1994). Common cladistic information and its consensus
          representation: reduced Adams and reduced cladistic consensus trees
          and profiles. Systematic Biology, 43(3), 343-368.

    Parameters
    ----------
    infilelist : list of str
        Paths of the file containing:

            * a hierarchical matrix (see 
              https://vrineau.github.io/AgattaDocs/Input%20files.html
              for informations about the format)
            * a single newick tree on each line. Example:

            (a,b,(c,d,(e,f)));
            ((c,e),(a,(b,(d,f))));

        The default is False (a selection window appears in this case).
    software_path : TYPE
        Path of the software designated with the "software" argument.
    software : str, optional
        Choose the software to use between 'paup', 'tnt', 'wqfm', 
        and 'wtree-qmc'.
        The default is 'paup'.
    taxa_replacement : str, optional
        Path of a table file containing two columns. The first column
        corresponds to the names of the terminals of the newick stored in
        infile, and the second column corresponds to their names the user wants
        to obtain at the end. All separators accepted. Example:
             AA Diceras
             AB Valletia
             AC Monopleura
        The default is False (no replacement).
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    weighting : str
        Weighting scheme to use between:

            * FW: Fractional weighting from Rineau et al. (2021)
            * FWNL: Fractional weighting from Nelson and Ladiges (1992)
            * UW (Uniform weighting from Nelson and Ladiges 1992)
            * MW: Minimal weighting from Wilkinson et al. (2004)
            * AW: Additive weighting : the weight of a triplet in additive
              weighting corresponds to the number of trees in which the
              triplet is present.
            * NW: No weighting (all triplets have a weight of 1).

        Nelson, G., & Ladiges, P. Y. (1992). Information content and fractional
        weight of three-item statements. Systematic biology, 41(4), 490-494.
        Rineau, V., Zaragüeta, R., & Bardin, J. (2021). Information content of
        trees: three-taxon statements, inference rules and dependency.
        Biological Journal of the Linnean Society, 133(4), 1152-1170.
        Wilkinson, M., Cotton, J. A., & Thorley, J. L. (2004). The information
        content of trees and their matrix representations.
        Systematic Biology, 53(6), 989-1001.

        The default is 'FW'.
    parallel : str or int
        Option for choosing if the analysis is made using multiprocessing or
        not. This argument can be:

            * "not" if the user does not wan to use multiprocessing.
            * "auto" for automatic detection of the number of cpu available.
            * any integer corresponding to the number of cpu the user wants to
              allow to the analysis.

    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. If set to none, a window appears to choose the location and name
        of the prefix for the nexus file. The default is 'agatta_out'.
    analysis : str
        Add a line in the nexus file to state the type of analysis between
        branch and bound ("bandb") or heuristic ("heuristic").
        The default is 'heuristic'.
    nrep : int, optional
        Add a line in the file to state the number of replicates in case
        of heuristic search. The default is 1000.
    biogeo_tab : str
        path of table with semicolon separators and two columns, taxa in
        the left column and corresponding areas in the right column.
        The default is False (no standardisation is made).
    chartest : bool, optional
        If True, write two files with the results of the character state test
        procedure. The default is False.
    ri : bool, optional
        If True, write a file with a global retention index and one retention
        index for each character. The default is False.
    consensus : bool or str, optional
        Computation of a consensus tree between:

            * strict: strict consensus
            * rcc: reduced cladistic consensus (Wilkinson 1994).
            * False: no cosnensus computed.

        The default is False.
    pdf_file : bool or str, optional
        Can be used to compute pdf files if the character state procedure.
        pdf_file argument is a path to a floder where Agatta will write one pdf
        file for each character state. The default is False (no pdf generated).
    dec_detail : bool
        If true, save a detailed output table in csv format named 
        prefix.triplet_table.csv of triplet weights per character in addition 
        to the .triplet file. 
    verbose : bool, optional
        Verbose mode if True. The default is True.

    Returns
    -------
    None.

    """

    print("Starting analysis")

    with open(prefix+".log", "w") as log_file:
        log_file.write("Agatta parameters\n\n")
        log_file.write("Current date and time : ")
        now = datetime.datetime.now()
        log_file.write(now.strftime("%Y-%m-%d %H:%M:%S"))
        log_file.write("\n")
        log_file.write("Input file paths: "+ ", ".join(infilelist)+"\n")
        log_file.write("Prefix: "+prefix+"\n")
        log_file.write("Taxa replacement file:"+str(taxa_replacement)+"\n")
        log_file.write("Method: "+method+"\n")

        if rosetta:
            log_file.write("Standardisation: yes\n")
            log_file.write("Table path for standardisation: "+rosetta+"\n")
        else:
            log_file.write("Standardisation: no\n")

        log_file.write("Software used: "+software+"\n")

        if software_path:
            log_file.write("Software path: "+software_path+"\n")
        else:
            log_file.write("Software path: no file path\n")

        log_file.write("Weighting triplets: "+weighting+"\n")
        log_file.write("Parallelisation: "+parallel+"\n")

        if software == "wqfm" or software == "wtree-qmc": # always heuristic
            log_file.write("Analysis: heuristic")
        elif analysis == "heuristic":
            log_file.write("Analysis: heuristic with {} replicates\n".format(
            str(nrep)))
        else:
            log_file.write("Analysis: " + analysis + "\n")

        log_file.write("Character states test: "+str(chartest)+"\n")

        if pdf_file:
            log_file.write("Character states test pdf results location: " +
                           str(pdf_file)+"\n")
        else:
            log_file.write("Character states test: no pdf required by user\n")

        log_file.write("Retention index: "+str(ri)+"\n")
        log_file.write("Consensus: "+str(consensus)+"\n")
        log_file.write("Verbose: "+str(verbose)+"\n")

    # if input is a matrix, use hmatrix to convert in a tree list
    character_dict = wrapper_character_extraction(infilelist, taxa_replacement,
                                                  prefix, verbose)

    # standardisation option
    if rosetta:

        #save character_dict
        character_dict = standardisation(character_dict,
                                         rosetta,
                                         prefix,
                                         verbose=verbose)

    # remove automatically repetitions if detected (user message printed)
    if rep_detector(character_dict):
        character_dict = del_replications_forest(character_dict,
                                                 method=method,
                                                 prefix=prefix,
                                                 verbose=verbose)

    prefix_path = os.path.join(os.getcwd(), prefix)

    # compute triplets and weights from the tree list and save a nex/tnt file
    convert(character_dict,
            "trees",
            prefix,
            parallel,
            weighting,
            analysis,
            taxa_replacement,
            nrep,
            True,  # logfile
            software,
            dec_detail,
            method,
            verbose)

    # three-item analysis using PAUP*, TNT, WQFM, or wTREE-QMC
    if software == "paup":
        prefix_end = prefix_path+".nex"
        search_pipeline(prefix_end,
                        software_path,  # except if software_path == False
                        "paup")

    elif software == "tnt":
        search_pipeline(prefix+".tnt",
                        software_path,
                        "tnt")

    elif software == "wqfm":
        search_pipeline(prefix+".wqfm",
                        software_path,
                        "wqfm", prefix)

    elif software == "wtree-qmc":
        search_pipeline(prefix+".wtqmc",
                        software_path,
                        "wtree-qmc", prefix)
        
    # extract trees from output file and delete root
    cladogram_dict = dict()
    i = 1

    with warnings.catch_warnings():  # rerooting func poorly tested
        warnings.filterwarnings("ignore", category=UserWarning)
        if software == "tnt":

            with open(prefix + ".tre", "r") as tntfile:
                tntfile2 = []
                newickstring = False
                for line in tntfile:
                    if "[&U] \n" in line or "[&R] \n" in line:
                        tntfile2.append(line[:-2])  #remove tnt \n
                        newickstring = True
                    elif newickstring:
                        tntfile2.append(line.replace(" ",""))
                        newickstring = False
                    else:
                        tntfile2.append(line)
            with open(prefix + ".tre", "w") as tntfile:
                for line in tntfile2:
                    tntfile.write(line)

            tstreelist = treeswift.read_tree_nexus(prefix + ".tre")

        elif (software == "paup" or software == "wqfm" or 
              software == "wtree-qmc"):
            try:
                tstreelist = treeswift.read_tree_newick(prefix + ".tre")
            except RuntimeError:
                print("ERROR: " + software + " failed to find a tree.")
                if prefix[0] == "_" and software == "paup":
                    print("PAUP* doesn't like prefix names starting with _")
                sys.exit(1)

        if not isinstance(tstreelist, list):
            tstreelist = [tstreelist]

        for tstree in tstreelist:
            if isinstance(tstree, dict):
                for idt, tst in tstree.items():
                    nodedict = tst.label_to_node(selection='leaves')
                    tst.reroot(nodedict["root"])
                    tst.suppress_unifurcations()
                    newickstring = tst.newick().replace("root","")

                    if (newickstring.startswith("[&R] ") or
                        newickstring.startswith("[&U] ")):
                        cladogram_dict[Tree(
                            newickstring.split(" ")[1])] = idt
                    else:
                        cladogram_dict[Tree(newickstring)] = idt

            else:
                nodedict = tstree.label_to_node(selection='leaves')
                try:
                    tstree.reroot(nodedict["root"])
                except KeyError:
                    print("ERROR: A problem occurred when rooting. "
                          + "The problem may be due to old output files left " 
                          + "in the working directory. Please try to delete"
                          + " all files in the folder other than input files "
                          + "and restart the analysis. ")
                    
                    sys.exit(1)
                tstree.suppress_unifurcations()
                newickstring = tstree.newick().replace("root","")

                if (newickstring.startswith("[&R] ") or
                    newickstring.startswith("[&U] ")):
                    cladogram_dict[Tree(newickstring.split(" ")[1])] = i
                else:
                    cladogram_dict[Tree(newickstring)] = i
                i += 1

    # branch length computation
    triplet_dict = triplet_extraction(prefix+'.triplet', 
                                      taxa_replacement_file=prefix+'.taxabloc')
    
    for tree in cladogram_dict.keys():
        for node in tree.traverse():
            if not node.is_leaf() and not node.is_root():
                
                node.dist = 0  # set support value to 0
                
                # for each node, compute support values
                for trip, FW in triplet_dict.items():
                    
                    in1, in2 = trip.in_taxa    
                    out, = trip.out_taxa    
                    
                    if (in1 in node.get_leaf_names() 
                        and in2 in node.get_leaf_names() 
                        and out not in node.get_leaf_names()):
                        
                        node.dist += FW

    # output file
    with open(prefix + ".tre", "w") as result_file:
        for tree in cladogram_dict.keys():
            result_file.write(tree.write(format=6) + "\n") #with support values

    if len(cladogram_dict.keys()) == 1:
        print("1 optimal tree found")
    else:
        print(str(len(cladogram_dict.keys())) + " optimal trees found")
        
    # consensus computation
    if consensus:

        # strict consensus
        if consensus == "strict":

            constrict(list(cladogram_dict.keys()), prefix)

        # reduced cladistic consensus (Wilkinson, 1994)
        elif consensus == "rcc":

            rcc(list(cladogram_dict.keys()), prefix)

    # character states test on the strict consensus
    if chartest:
        character_states_test({constrict(list(cladogram_dict.keys()),
                                         silent=True): 1},
                              character_dict,
                              prefix,
                              pdf_file)

    # retention index
    if ri:
        RI(cladogram_dict,
           character_dict,
           weighting=weighting,
           prefix=prefix)

    print('The analysis ended successfully')
