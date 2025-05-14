# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from .ini import wrapper_character_extraction
from .ini import taxa_to_numbers
from .ini import infotree_checker
from fractions import Fraction
from itertools import combinations
from itertools import product
from multiprocessing import Pool
from multiprocessing import cpu_count
from collections import defaultdict
from functools import partial
from tqdm import tqdm
import os
import sys
import csv
import time
import uuid
import shutil
import pickle
import pandas as pd
from ete3 import Tree


class triplet():

    __slots__ = ("in_taxa", "out_taxa", "FW")

    def __init__(self, in_taxa=set(), out_taxa=set(), FW=Fraction()):

        self.in_taxa = set(in_taxa)  # taxa in the informative node
        self.out_taxa = set(out_taxa)  # sister group of the informative node
        self.FW = FW  # Triplet weight

    def __eq__(self, other_triplet):

        return (set(self.in_taxa) == set(other_triplet.in_taxa)) and (
            self.out_taxa == other_triplet.out_taxa)

    def __call__(self, in_taxa=set(), out_taxa=set(), FW=Fraction(),
                 parent_triplets=()):

        return self

    def __hash__(self):

        return hash(list(self.in_taxa)[0]) ^ hash(list(self.in_taxa)[1]) ^ \
            hash(list(self.out_taxa)[0])

    def __del__(self):

        del self

    def __repr__(self):

        return "({},({},{}))".format(str(list(self.out_taxa)[0]),
                                    str(sorted(list(self.in_taxa))[0]),
                                    str(sorted(list(self.in_taxa))[1]))


def del_replications(treerep, method="TMS", verbose=False):
    """
    Remove all repeated leaves from a single rooted tree according to the
    free-paralogy subtree analysis. Two algorithms are currently implemented in
    Agatta:

        * The original algorithm of Nelson and Ladiges (1996) designed in the
          paradigm of cladistic biogeography: Free-paralogy Subtree analysis.

        * The algorithm of Rineau et al. (2022) designed to construct
          subtrees without repeated leaves while minimising the loss of
          information in terms of triplets: Triplet Maximisation Subtrees.
          This is the default mode.

          Nelson, G. J., & Ladiges, P. Y. (1996). Paralogy in cladistic
          biogeography and analysis of paralogy-free subtrees.
          American Museum novitates 3167.

          Rineau, V., Moncel, M-H., & Zeitoun, V. (2022). Revealing
          evolutionary patterns behind homogeneity: the case of the Paleolithic
          assemblages from Notarchirico (Southern Italy).

    Parameters
    ----------
    treerep : ete3.Tree
        One rooted tree.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    verbose : bool, optional
        Verbose mode if true. The default is False.

    Returns
    -------
    list
        A list of subtrees without repeated leaves.

    """

    def del_replications_node(treerep, method="TMS", verbose=False):
        """
            Function to remove a single repeated leaf.
            Used in del_replications.
        """

        sys.setrecursionlimit(10**7)  # can hit the recursion limit

        paralogs = dict()
        cardl = list(Tree.get_leaf_names(treerep))
        cards = set(Tree.get_leaf_names(treerep))
        cardr = [x for x in cards if cardl.count(x) > 1]  # list of rep. leaves
        treelistrep = []
        treelistnorep = []
        treelistverif = []

        if cardr:  # if repetitions
            nodelist = sorted([node for node in treerep.traverse(
                strategy="postorder") if not node.is_leaf()],
                key=lambda n: len(n))

            for e in nodelist:  # iterate from the smallest to the biggest node
                leafset = set(Tree.get_leaf_names(e))

                # analysis of the smallest node with repetitions
                if len(leafset) != len(Tree.get_leaf_names(e)):  # if rep.

                    # Identification if rep. is: apical, ortholog, paralog
                    # apical: node including only leaves
                    # ortholog: asymmetric node
                    # paralog: symmetric node
                    i = 0  # number of internal nodes connected to e
                    e.add_feature("FP", "main")
                    for child_node in e.get_children():
                        if not child_node.is_leaf():
                            i += 1
                            for cnode in child_node.traverse(
                                    strategy="postorder"):
                                cnode.add_feature("FP", i)  # label nodes save

                    for node in treerep.traverse():
                        try:
                            if node.FP:
                                pass
                        except AttributeError:
                            node.add_feature("FP", "out")

                    # special case: repeated leaf branched to a symmetric node
                    paralog_but_leaf = False

                    if i > 1:
                        special_tree = treerep.copy(method='cpickle')
                        for child in special_tree.search_nodes(
                                FP="main")[0].get_children():
                            if child.is_leaf:
                                if child.name in cardr:
                                    paralog_but_leaf = True
                                    if verbose:
                                        print("ortholog repetition (on "
                                              + "symmetric node):{}".format(
                                                  child.name))
                                    child.delete()

                    # special case: leaves on paralog
                    if paralog_but_leaf:
                        paralogs["main"] = special_tree

                    # apical
                    elif i == 0:
                        paralogs["main"] = treerep.copy(method='cpickle')
                        for l in set(cardr) & leafset:  # for each leaf
                            leaf1 = False
                            for delnode in paralogs["main"].search_nodes(
                                    FP="main")[0].get_leaves_by_name(name=l):
                                if not leaf1:
                                    leaf1 = True  # except the first
                                else:
                                    delnode.delete()  # delete leaves
                                    if verbose:
                                        print("apical repetition:{}".format(
                                            delnode.name))

                    # ortholog
                    elif i == 1:
                        paralogs["main"] = treerep.copy(method='cpickle')
                        kl = []

                        for l in set(cardr) & leafset:  # for each leaf

                            leaf1 = False

                            for sorted_nodes in sorted(
                                    [node for node in paralogs[
                                        "main"].search_nodes(
                                            FP="main")[0].traverse(
                                                strategy="postorder") if not
                                        node.is_leaf()], key=lambda n: len(n)):

                                for childnode in sorted_nodes.get_children():
                                    if childnode.is_leaf():
                                        if childnode.name == l:  # for each
                                            if not leaf1:  # rep. leaf
                                                leaf1 = True    # except first
                                            elif leaf1:
                                                kl.append(childnode)

                        for delnode in kl:
                            delnode.delete()  # delete leaves
                            if verbose:
                                print("ortholog repetition:{}".format(
                                    delnode.name))

                    # paralog
                    else:

                        j = 1

                        # Rineau's algorithm - Triplet Maximisation Subtrees
                        if method == "TMS":
                            paralogs["main"] = treerep.copy(method='cpickle')

                            for delnode in paralogs["main"].search_nodes(
                                   FP="main")[0].iter_descendants("postorder"):
                                if not delnode.is_leaf():
                                    delnode.delete()  # del internal node

                            while j != i+1:
                                paralogs[j] = treerep.copy(method='cpickle')

                                for delnode in paralogs[j].traverse():
                                    if (not delnode.is_leaf() and not
                                    delnode.FP == j):
                                        delnode.delete()  # del internal node

                                j += 1

                        # Nelson's algorithm - Free paralogy subtree analysis
                        elif method == "FPS":

                            while j != i+1:  # for each node brch to paralog

                                # subtree for each node branched to the paralog
                                paralogs[j] = treerep.copy(method='cpickle')

                                for delnode in paralogs[j].search_nodes(
                                        FP="main")[0].get_children():
                                    if not delnode.FP == j:
                                        delnode.detach()

                                paralogs[j].search_nodes(FP="main")[0].delete()
                                j += 1

                        if verbose:
                            print("paralog repetition. " + str(len(paralogs)) +
                                  " subtrees " +
                                  "generated from the main tree.")

                    break  # function only for one instance of repetition

            treelistverif = [p for p in paralogs.values()]

        else:  # if there is no repetition
            treelistnorep.append(treerep)

        # sort between trees with and without repetitions
        for t in treelistverif:
            if len(set(Tree.get_leaf_names(t))) > 2:
                if len(list(Tree.get_leaf_names(t))) != len(set(
                        Tree.get_leaf_names(t))):  # detect repetitions
                    treelistrep.append(t)
                else:
                    treelistnorep.append(t)

        # del features (main) in all trees for recursion
        treelistrepclearcopy = []
        for t in treelistrep:
            for node in t.traverse():
                node.del_feature('FP')
            treelistrepclearcopy.append(t.copy(method='cpickle'))

        return treelistrepclearcopy, treelistnorep

    listundone, listdone = del_replications_node(treerep, method)

    if listundone:  # if repetitions
        while listundone:  # while until listundone contains trees to analyse
            listundone1 = []
            listdone1 = []

            for treerep2 in listundone:
                LU, LD = del_replications_node(treerep2, method)
                listundone1 += LU
                listdone1 += LD

            # all trees analysed
            listundone = listundone1.copy()  # replace trees with rep. by list
            listdone += listdone1  # accumulation of finished trees in listdone

    # detect non-informative trees
    treelist = []
    for infotree in listdone:

        # detect and delete empty internal node (one single descendant)
        for node in infotree.traverse():
            if not node.is_leaf() and len(node.get_children()) == 1:
                node.get_children()[0].delete()

        # retain tree only if more than one internal node
        cdict, noninfo = infotree_checker({infotree:1}, verbose=False)
        if cdict.keys():
            infotree = list(cdict.keys())[0]
            infotree.ladderize()
            treelist.append(infotree)

    if verbose:  # verbose mode

        print("Output trees:")
        if treelist:
            for l in treelist:
                print(l.write(format=9))
        else:
            print("No informative tree")
        print("")

    return treelist


def del_replications_forest(character_dict, method="TMS",
                            prefix="agatta_del_replications", verbose=False):
    """
    Remove all repeated leaves from trees according to the
    free-paralogy subtree analysis. Two algorithms are currently implemented in
    Agatta:

        * The original algorithm of Nelson and Ladiges (1996) designed in the
          paradigm of cladistic biogeography.

        * The algorithm of Rineau et al. (2022) designed to construct
          subtrees without repeated leaves while minimising the loss of
          information in terms of triplets. This is the default mode.

          Nelson, G. J., & Ladiges, P. Y. (1996). Paralogy in cladistic
          biogeography and analysis of paralogy-free subtrees.
          American Museum novitates 3167.

          Rineau, V., Moncel, M-H., & Zeitoun, V. (2022). Revealing
          evolutionary patterns behind homogeneity: the case of the Paleolithic
          assemblages from Notarchirico (Southern Italy).

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    prefix : str, optional
        Prefix of the saving file. The default is False.
    verbose : bool, optional
        Verbose mode if true. The default is False.

    Returns
    -------
    tree_dict
        A dictionary of resulting subtrees in keys and their id as values.
        Id is a string. If a tree 2 results in a single subtree, the id of the
        subtree will be '[2]', otherwise if several subtrees are generated
        their respective id will be '[2.1]', '[2.2]', etc.

    """
    print("Managing polymorphism")
    tree_dict = dict()
    output_trees = 0


    if verbose:
        loopchar = character_dict.items()
    else:
        loopchar = tqdm(character_dict.items())
        
    if prefix:
        with open(prefix+".poly", "w") as logfile:
            logfile.write("Polymorphism management: {}\n".format(method))

    for treed, treeid in loopchar:

        if verbose:
            print("[" + str(treeid) + "] " + treed.write(format=9))
        
        if prefix:
            with open(prefix+".poly", "a") as logfile:
                logfile.write("[" + str(treeid) + "] " 
                              + treed.write(format=9) + "\n")

        treelist = del_replications(treed, method, verbose)

        if treelist:
            
            # if single tree obtained
            if len(treelist) == 1:
                tree_dict[treelist[0]] = treeid
                
                if prefix:
                    if len(treed) == len(treelist[0]):
                        with open(prefix+".poly", "a") as logfile:
                            logfile.write("[" + str(treeid) + "] " +
                                          "no repeated leaf\n")
                    else:
                        with open(prefix+".poly", "a") as logfile:
                            logfile.write("[" + str(treeid) + "] " +
                                          treelist[0].write(format=9) + "\n")
                        
            # if character cut into several free-paralogy subtrees
            else:
                
                # function to sort trees with their minimal state number 
                def sorttreelist(character):
                    """
                    function to sort trees with their min state nb
                    return empty list otherwise
                    """
                    statenb = []
                    
                    try: 
                        for n in character.traverse(strategy="preorder"):
                            if (n.is_leaf() == False and 
                                n.is_root() == False):
                                statenb.append(int(n.charstate_name))
                        
                        return min(statenb) 
                    
                    except AttributeError:
                        return 1000
                          
                #reorder list with state numbers if state numbers
                treelist.sort(key=sorttreelist)

                i = 1
                output_trees += len(treelist)
                for treel in treelist:
                    tree_dict[treel] = str(treeid) + "." + str(i)


                    if prefix:
                        with open(prefix+".poly", "a") as logfile:
                            logfile.write("[" + str(treeid) + "." + str(i) + 
                                          "] " + treel.write(format=9) + "\n")
                    i += 1
        else:
            if prefix:
                with open(prefix+".poly", "a") as logfile:
                    logfile.write("[" + str(treeid) + 
                                  "] no informative tree\n")

        if prefix:
            with open(prefix+".poly", "a") as logfile:
                logfile.write("\n")

    td, noninfo = infotree_checker(tree_dict)
    
    # message when non-informative characters
    for n in noninfo:
        print("Character {} non-informative after " 
              +"polymorphism management".format(n))

    if output_trees != 0:
        print("Polymorphism removed: "
              "{} polymorphism free informative characters "
              "computed from {} character trees".format(
                  str(output_trees), str(len(character_dict))))
    else:
        print("Polymorphism removed. No polymorphism free "
              "informative subtree remaining.")
        sys.exit(1)

    return tree_dict


def rep_detector(character_dict):
    """
    Function to detect if repeated leaves are present in a dictionary of trees.

    """

    for t in character_dict.keys():
        cardl = list(Tree.get_leaf_names(t))
        cards = set(Tree.get_leaf_names(t))
        if len(cardl) != len(cards):
            return True

    return False


def triplet_extraction(infile, taxa_replacement_file=False):
    """
    Extract a dictionnary of triplets (keys) and weights (values) from a 
    prefix.triplet file generated by main_tripdec. Can use a taxa bloc 
    generated by main_tripdec for taxa name replacement.

    Parameters
    ----------
    infile : str
        Path of file containing triplets and weights.
    taxa_replacement_file : str, optional
        Path of a table file containing two columns. The first column
        corresponds to the names of the terminals of the newick stored in
        infile which must be integers, and the second column corresponds to
        their names the user wants to obtain at the end.
        All separators accepted. Example:
             1 Diceras
             2 Valletia
             3 Monopleura
        The default is False (no replacement).
    Returns
    -------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values).

    """

    print("Loading triplet set")

    triplet_dict = {}

    infile = os.path.expanduser(infile)

    if not os.path.isfile(infile):
        print("ERROR: The file '" + infile + "' does not exist." +
                       "\nOperation aborted.")
        sys.exit(1)

    if taxa_replacement_file:
        taxa_dict = {}
        # detect separator in taxa_replacement_file
        with open(taxa_replacement_file, "r") as taxa_table1:
            try:
                dialect = csv.Sniffer().sniff(taxa_table1.read())
            except:
                print("ERROR: Could not determine separator in the "
                                  + "file '" + taxa_replacement_file
                                  + "'. The table is probably broken."
                                  + "\nOperation aborted.")
                sys.exit(1)

        # build taxanames to id converter
        with open(taxa_replacement_file, "r") as taxa_table1:
            data = csv.reader(taxa_table1, delimiter=dialect.delimiter)
            taxa_table2 = list(data)

            rowlist = [[i for i in row if i] for row in taxa_table2]

            for idtax, nametax in rowlist:
                taxa_dict[int(idtax)] = nametax

    # read triplet file and build triplet dictionary
    with open(infile, "r") as file_tree:
        for line in file_tree:
            if line.strip():
                tripletstr = line.split("    ")
                newstr = ''.join(
                    (ch if ch in '0123456789.-e' else ' ')
                    for ch in tripletstr[0])
                taxaint = [int(i) for i in newstr.split()]

                if taxa_replacement_file:
                    for i in [0,1,2]:
                        try:
                            taxa_dict[taxaint[i]]
                        except KeyError:
                            print("ERROR: The name '"
                                           + str(taxaint[i])
                                           + "' does not exists in the table "
                                           + "file '"
                                           + taxa_replacement_file
                                           + "'\nOperation aborted.")
                            sys.exit(1)

                    convert_trip = triplet({taxa_dict[taxaint[1]],
                                            taxa_dict[taxaint[2]]},
                                           {taxa_dict[taxaint[0]]})
                    triplet_dict[convert_trip] = Fraction(tripletstr[1])

                else:
                    trip = triplet({taxaint[1], taxaint[2]}, {taxaint[0]})
                    triplet_dict[trip] = Fraction(tripletstr[1])

    print("{} triplets loaded".format(str(len(triplet_dict))))

    return triplet_dict


def picklemerge(namelist):
    """
    Take as input a list of pickle file names, extract the triplet dictionaries
    and merge them.

    Parameters
    ----------
    namelist : list
        List of file names.

    Returns
    -------
    superdict : dict
        Dictionary of triplets (keys) and weights (values).

    """

    # load first pickle (first component triplet set)
    with open(namelist[0], 'rb') as pickle_file:
        superdict = defaultdict(int, pickle.load(pickle_file))

    # load each pickle file and add triplets to dict / cumulative weights
    for i in range(1, len(namelist)):
        with open(namelist[i], 'rb') as pickle_file:
            tripletdict_temp = pickle.load(pickle_file)
            for trip, FW in tripletdict_temp.items():
                superdict[trip] += FW

    return superdict


def tripdecFW(triplet_output=dict(), total_taxaset=False, character=Tree()):
    """
    Decompose a single tree into its triplets and compute their weights using
    the fractional weighting following Rineau et al. (2021).

    Rineau, V., Zaragüeta, R., & Bardin, J. (2021). Information content of
    trees: three-taxon statements, inference rules and dependency.
    Biological Journal of the Linnean Society, 133(4), 1152-1170.

    Parameters
    ----------
    triplet_output : dict, optional
        Used for recursion. The default is dict().
    total_taxaset : TYPE, optional
        Used for recursion. The default is False.
    character : ete3.Tree, optional
        A tree to decompose into triplets. The default is Tree().

    Returns
    -------
    dict
        Dictionary of triplets (keys) and weights (values).

    """

    if not total_taxaset:
        total_taxaset = set(character.get_leaf_names())

    # Compute triplet dictionary
    tree_tripdic = defaultdict(int)  # dict of triplets from tree
    nodaltripletdict = defaultdict(int)  # dict of triplets from component
    children_generator = (child_node for child_node in character.get_children(
        ) if not child_node.is_leaf())  # compute descendants info. nodes

    # Compute triplets from descendant nodes
    for child_node in children_generator:
        newtripletdict = tripdecFW(defaultdict(int),
                                   total_taxaset,
                                   character=child_node)

        for trip, FW in newtripletdict.items():
            tree_tripdic[trip] += FW

    # Compute triplets if not is not root
    if character.is_root():  # if node is root, end of analysis

        return tree_tripdic  # end

    else:
        taxa_in = set(character.get_leaf_names())  # taxa inside the node
        taxa_out = total_taxaset - taxa_in  # taxa outside the node
        totalnodeweight = len(taxa_out)*(len(taxa_in)-1)

        for taxa_in1, taxa_in2 in combinations(taxa_in, 2):  # all leaf couples
            for singleout in taxa_out:  # all sister taxa
                if not triplet({taxa_in1, taxa_in2},
                               {singleout}) in tree_tripdic:
                    nodaltripletdict[triplet({taxa_in1,
                                              taxa_in2},
                                             {singleout})] = 0

        for trip in nodaltripletdict.keys():
            tree_tripdic[trip] += Fraction(totalnodeweight,
                                           len(nodaltripletdict))

        return tree_tripdic


def tripdec(weighting, character):
    """
    Decompose a single tree into its triplets and compute their weights using
    a weighting scheme between:

        * Fractional weighting from Nelson and Ladiges (1992)
        * Uniform weighting from Nelson and Ladiges (1992)
        * Minimal weighting from Wilkinson et al. (2004)
        * Additive weighting : the weight of a triplet in additive weighting
        corresponds to the number of trees in which the triplet is present.
        * No weighting (all triplets have a weight of 1).

        Nelson, G., & Ladiges, P. Y. (1992). Information content and fractional
        weight of three-item statements. Systematic biology, 41(4), 490-494.

        Wilkinson, M., Cotton, J. A., & Thorley, J. L. (2004). The information
        content of trees and their matrix representations.
        Systematic Biology, 53(6), 989-1001.

    Parameters
    ----------
    weighting : str
        Weighting scheme to use between: * FWNL (Fractional weighting from
                                               Nelson and Ladiges),
                                         * UW (Uniform weighting from Nelson
                                               and Ladiges)
                                         * MW (Minimal Weighting),
                                         * AW (Additive Weighting),
                                         * NW (No Weighting).
    character : ete3.Tree
        A tree to decompose into triplets. The default is Tree().

    Returns
    -------
    dict
        Dictionary of triplets (keys) and weights (values).

    """

    triplet_output = dict()
    cardinal_character1 = set(Tree.get_leaf_names(character))  # leaf names
    internal_nodes1 = (node for node in character.traverse(
        strategy="preorder") if not node.is_leaf() and not node.is_root())
    tree_tripdic = defaultdict(int)
    card = len(cardinal_character1)
    MW_node = Fraction(6, card * (card - 1))

    for node in internal_nodes1:
        taxa_in = set(Tree.get_leaf_names(node))  # taxa inside the node
        taxa_out = cardinal_character1 - taxa_in  # taxa outside the node
        FW_node = Fraction(2, len(taxa_in))

        # build triplets from the node
        for taxa_in1, taxa_in2 in combinations(taxa_in, 2):
            for singleout in taxa_out:
                if weighting == "FWNL":
                    tree_tripdic[triplet({taxa_in1, taxa_in2},
                                         {singleout})] += FW_node

                elif weighting == "UW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2},
                                         {singleout})] += 1

                if weighting == "MW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2},
                                         {singleout})] = MW_node

                elif weighting == "AW" or weighting == "NW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2},
                                         {singleout})] = 1

    if weighting == "NW":  # if no weighting

        triplet_output.update(tree_tripdic)

        return triplet_output

    else:
        for trip in tree_tripdic.keys():
            if trip in triplet_output:  # add new triplets
                triplet_output[trip] += tree_tripdic[trip]
            else:  # add weight if triplet already exists
                triplet_output[trip] = tree_tripdic[trip]

        return triplet_output


def tripdec_allweights(weighting, character):
    """
    Function to be used in parallel_tripdec only for multiprocessing
    computation.
    Decompose a single tree into its triplets and compute their weights
    using a weighting scheme between:

        * Fractional weighting from Rineau et al. (2021)
        * Fractional weighting from Nelson and Ladiges (1992)
        * Uniform weighting from Nelson and Ladiges (1992)
        * Minimal weighting from Wilkinson et al. (2004)
        * Additive weighting : the weight of a triplet in additive weighting
        corresponds to the number of trees in which the triplet is present.
        * No weighting (all triplets have a weight of 1).

        Nelson, G., & Ladiges, P. Y. (1992). Information content and fractional
        weight of three-item statements. Systematic biology, 41(4), 490-494.

        Rineau, V., Zaragüeta, R., & Bardin, J. (2021). Information content of
        trees: three-taxon statements, inference rules and dependency.
        Biological Journal of the Linnean Society, 133(4), 1152-1170.

        Wilkinson, M., Cotton, J. A., & Thorley, J. L. (2004). The information
        content of trees and their matrix representations.
        Systematic Biology, 53(6), 989-1001.

    Parameters
    ----------
    weighting : str
        Weighting scheme to use between: * FW (Fractional weighting from
                                               Rineau et al. 2021)
                                         * FWNL (Fractional weighting from
                                               Nelson and Ladiges),
                                         * UW (Uniform weighting from
                                               Nelson and Ladiges (1992)
                                         * MW (Minimal Weighting),
                                         * AW (Additive Weighting),
                                         * NW (No Weighting).
    character : ete3.Tree
        A tree to decompose into triplets. The default is Tree().

    Returns
    -------
    pickle_dict : dict
        Dictionary of triplets (keys) and weights (values) used for parallel
        computation in parallel_tripdec.
    """

    if not weighting in ("FW", "FWNL", "UW", "MW", "AW", "NW"):
        print("ERROR: " + weighting +
                       " is not a correct weighting method.\n" +
                       "Allowed weightings: FW, FWNL, MW, UW, AW, NW" +
                       "\nOperation aborted.")
        sys.exit(1)

    pickle_dict = defaultdict(int)
    tree_id = str(uuid.uuid4())

    cardinal_character1 = set(Tree.get_leaf_names(character))  # leaf names
    internal_nodes1 = (node for node in character.traverse(
        strategy="postorder") if not node.is_leaf() and not node.is_root())
    card = len(cardinal_character1)
    MW_node = Fraction(6, card * (card - 1))

    for node in internal_nodes1:
        pickle_name = str(uuid.uuid4()) + "_agatta.pickle"
        nodaltripletdict = dict()
        taxa_in = set(Tree.get_leaf_names(node))  # taxa inside the node
        taxa_out = cardinal_character1 - taxa_in  # taxa outside the node
        FWNL_node = Fraction(2, len(taxa_in))
        tempw = sum((len(t1)*len(t2) for t1, t2 in combinations((
            child_node.get_leaf_names() for child_node in node.get_children(
                )), 2)))
        FW_node = Fraction(len(taxa_out)*(len(taxa_in)-1), len(taxa_out)*tempw)

        for singleout in taxa_out:

            # generate all exclusive triplets from component
            if weighting in ("FW", "MW", "AW", "NW"):

                for taxalist1, taxalist2 in combinations((
                        child_node.get_leaf_names(
                        ) for child_node in node.get_children()), 2):
                    for taxa_in1, taxa_in2 in product(taxalist1, taxalist2):

                        if weighting == "FW":
                            nodaltripletdict[triplet(
                                in_taxa={taxa_in1, taxa_in2},
                                out_taxa={singleout})] = FW_node

                        elif weighting == "MW":
                            nodaltripletdict[triplet(
                                in_taxa={taxa_in1, taxa_in2},
                                out_taxa={singleout})] = MW_node

                        else:  # if additive weighting or no weighting
                            nodaltripletdict[triplet(
                                in_taxa={taxa_in1, taxa_in2},
                                out_taxa={singleout})] = 1

            # generate all triplets from component
            elif weighting in ("FWNL", "UW"):

                for taxa_in1, taxa_in2 in combinations(taxa_in, 2):

                    if weighting == "FWNL":
                        nodaltripletdict[triplet(
                            in_taxa={taxa_in1, taxa_in2},
                            out_taxa={singleout})] = FWNL_node

                    else:  # UW
                        nodaltripletdict[triplet(
                            in_taxa={taxa_in1, taxa_in2},
                            out_taxa={singleout})] = 1

        with open(pickle_name, 'wb') as pickle_file:
            pickle.dump(nodaltripletdict, pickle_file,
                        protocol=pickle.HIGHEST_PROTOCOL)

        pickle_dict[pickle_name] = tree_id

    return pickle_dict


def standard_tripdec(character_dict, weighting, dec_detail=False, 
                     prefix=False, verbose=True):
    """
    Compute triplets and weights from several rooted trees.
    Each tree is decomposed one after the other.

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    weighting : str
        Weighting scheme to use between: * FW (Fractional weighting from
                                               Rineau et al. 2021)
                                         * FWNL (Fractional weighting from
                                               Nelson and Ladiges),
                                         * UW (Uniform weighting from
                                               Nelson and Ladiges 1992)
                                         * MW (Minimal Weighting),
                                         * AW (Additive Weighting),
                                         * NW (No Weighting).
    dec_detail : bool, optional
        If true, save a detailed table in csv format named 
        prefix.triplet_table.csv of triplet weights per character in addition 
        to the .triplet file. The default is False.
    prefix : str, optional
        Prefix of the saving file. The default is False.
    verbose : bool, optional
        Activate verbose mode if true. The default is True.

    Returns
    -------
    triplet_output : dict
        Dictionary of triplets (keys) and weights (values).

    """

    triplet_output = dict()
    triplet_output_per_char = dict()

    if weighting == "FW":

        for treedec, charnb in character_dict.items():
            triplet_output2 = tripdecFW(triplet_output,
                                        False, character=treedec)

            if dec_detail:
                triplet_output_per_char[charnb] = triplet_output2

            for trip, FW in triplet_output2.items():
                if trip in triplet_output:  # new triplets
                    triplet_output[trip] += FW
                else:  # sum triplet weights
                    triplet_output[trip] = FW

    elif weighting in ("FWNL", "MW", "AW", "UW", "NW"):
        for treedec, charnb in character_dict.items():
            triplet_output2 = tripdec(weighting, treedec)
            
            if dec_detail:
                triplet_output_per_char[charnb] = triplet_output2

            for trip, FW in triplet_output2.items():
                if trip in triplet_output and not weighting == "NW":
                    triplet_output[trip] += FW
                else:  # sum triplet weights
                    triplet_output[trip] = FW

    else:
        print("ERROR: " + weighting +
                       " is not a correct weighting method.\n" +
                       "Allowed weightings: FW, FWNL, MW, UW, AW, NW" +
                       "\nOperation aborted.")
        sys.exit(1)

    if dec_detail:
        
        # create detailed character table
        chartable = pd.DataFrame(0.000,
                         index=[str(trip) for trip in triplet_output.keys()], 
                         columns=triplet_output_per_char.keys())
        
        # fill table
        for charnb, triplet_output2 in triplet_output_per_char.items():
            for trip, FW in triplet_output2.items():
                
                chartable.loc[str(trip), charnb] += FW
                
        chartable.loc['Total'] = chartable.sum(axis=0)
        chartable.loc[:,'Total'] = chartable.sum(axis=1)
        chartable = chartable.round(decimals = 3)
        
        #save table
        chartable.to_csv(prefix + '.triplet_table.csv')
        print('Detailed triplet/character table computed.')
                    
    if prefix:
        with open(prefix+".triplet", "w") as tdfile:
            taxa_dict, code_dict, taxa_convers = taxa_to_numbers(
                                                                character_dict)
            if weighting in ("FW", "FWNL", "MW"):
                for trip, FW in triplet_output.items():

                    in_taxa = list(trip.in_taxa)
                    tdfile.write("(" + str(taxa_dict[list(trip.out_taxa)[0]])
                                 + ",("
                                 + str(taxa_dict[in_taxa[0]]) + ","
                                 + str(taxa_dict[in_taxa[1]]) + ")):    "
                                 + str(FW) + "    "
                                 + str(round(float(FW), 4))
                                 + "\n")

            else:
                 for trip, FW in triplet_output.items():

                    in_taxa = list(trip.in_taxa)
                    tdfile.write("(" + list(trip.out_taxa)[0] + ",("
                                 + in_taxa[0] + "," + in_taxa[1] + "));    "
                                 + str(FW) + "\n")

        with open(prefix+".taxabloc", "w") as taxa_bloc_file:
            for taxa, code in taxa_dict.items():
                taxa_bloc_file.write("{}    {}\n".format(str(code), taxa))

    return triplet_output


def parallel_tripdec(character_dict, weighting, prefix=False, ncpu="auto"):
    """
    Compute triplets and weights from several rooted trees in parallel.

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    weighting : str
        Weighting scheme to use between: * FW (Fractional weighting from
                                               Rineau et al. 2021)
                                         * FWNL (Fractional weighting from
                                               Nelson and Ladiges),
                                         * UW (Uniform weighting from
                                               Nelson and Ladiges 1992),
                                         * MW (Minimal Weighting),
                                         * AW (Additive Weighting),
                                         * NW (No Weighting).
    prefix : str, optional
        Prefix of the saving file. The default is False.
    ncpu : str, optional
        Number of cpu to be used. The default is "auto".

    Returns
    -------
    triplet_output : dict
        Dictionary of triplets (keys) and weights (values).

    """

    def split(a, n):
        """
            Sort in n lists a single list a of pickle names
        """
        k, m = divmod(len(a), n)
        return list(a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

    # number of cpu
    if ncpu == "auto":
        ncpu = cpu_count()
    else:
        ncpu = int(ncpu)

    tmp = "tmp_Agatta_picKles_" + str(uuid.uuid4())  # folder for pickles

    try:
        os.mkdir(tmp)
    except FileExistsError:
        shutil.rmtree(os.getcwd() + '/' + tmp, ignore_errors=True)
        os.mkdir(tmp)
    finally:
        os.chdir(tmp)

    # computing triplets for each tree in parallel
    treetuple = tuple(character_dict)  # tuple dict for multiprocessing

    # creating dictionary
    with Pool(ncpu) as pool:  # multiprocessing activation
        func = partial(tripdec_allweights, weighting)
        pickle_double = pool.map(func, treetuple)  # list of triplet dict

    # pickle dictionaries extraction and merging
    picknum = len([item for sublist in [tuple(dicname.keys(
        )) for dicname in pickle_double] for item in sublist])

    if ncpu > picknum:  # if more cpu than pickles, use minimal number of cpu
        ncpu = picknum
        pickletuple = split([item for sublist in [tuple(dicname.keys(
            )) for dicname in pickle_double] for item in sublist], ncpu)
    else:  # else, use the maximum number of cpu allowed
        pickletuple = split([item for sublist in [tuple(dicname.keys(
            )) for dicname in pickle_double] for item in sublist], ncpu)

    with Pool(ncpu) as pool:  # multiprocessing
        dictuple = pool.map(picklemerge, pickletuple)  # list of triplet dicts

    if weighting == "NW":  # unique weight of 1
        shared_dict = {t: 1 for t in dictuple[0]}
        for d in dictuple[1:len(dictuple)]:
            for trip in d:
                shared_dict[trip] = 1

    else:  # sum weights
        shared_dict = dictuple[0]
        for d in dictuple[1:len(dictuple)]:
            for trip, FW in d.items():
                shared_dict[trip] += FW

    # pickles deletion
    os.chdir('..')
    shutil.rmtree(os.getcwd() + '/' + tmp, ignore_errors=True)

    if prefix:
        with open(prefix+".triplet", "w") as tdfile:
            taxa_dict, code_dict, taxa_convers = taxa_to_numbers(
                                                                character_dict)
            if weighting in ("FW", "FWNL", "MW"):
                for trip, FW in shared_dict.items():

                    in_taxa = list(trip.in_taxa)
                    tdfile.write("(" + str(taxa_dict[list(trip.out_taxa)[0]])
                                 + ",("
                                 + str(taxa_dict[in_taxa[0]]) + ","
                                 + str(taxa_dict[in_taxa[1]]) + ")):    "
                                 + str(FW) + "    "
                                 + str(round(float(FW), 4))
                                 + "\n")

            else:
                 for trip, FW in shared_dict.items():

                    in_taxa = list(trip.in_taxa)
                    tdfile.write("(" + list(trip.out_taxa)[0] + ",("
                                 + in_taxa[0] + "," + in_taxa[1] + "));    "
                                 + str(FW) + "\n")

        with open(prefix+".taxabloc", "w") as taxa_bloc_file:
            for taxa, code in taxa_dict.items():
                taxa_bloc_file.write("{}    {}\n".format(str(code), taxa))

    return shared_dict


def main_tripdec(input_item, prefix, taxa_replacement, weighting, parallel, 
                 dec_detail, method, verbose):
    """
    Main function of Agatta for decomposing trees into triplets and compute
    triplet weights using multiprocessing or not. The input trees can be a
    dictionary containing ete3 trees as keys or the path(s) to input file(s)
    (https://vrineau.github.io/AgattaDocs/Input%20files.html).
    Decompose trees into their triplets and compute weights using a weighting
    scheme between:

        * Fractional weighting from Rineau et al. (2021)
        * Fractional weighting from Nelson and Ladiges (1992)
        * Uniform weighting from Nelson and Ladiges (1992)
        * Minimal weighting from Wilkinson et al. (2004)
        * Additive weighting : the weight of a triplet in additive weighting
        corresponds to the number of trees in which the triplet is present.
        * No weighting (all triplets have a weight of 1).

        Nelson, G., & Ladiges, P. Y. (1992). Information content and fractional
        weight of three-item statements. Systematic biology, 41(4), 490-494.

        Rineau, V., Zaragüeta, R., & Bardin, J. (2021). Information content of
        trees: three-taxon statements, inference rules and dependency.
        Biological Journal of the Linnean Society, 133(4), 1152-1170.

        Wilkinson, M., Cotton, J. A., & Thorley, J. L. (2004). The information
        content of trees and their matrix representations.
        Systematic Biology, 53(6), 989-1001.

    Parameters
    ----------
    input_item : string or dictionary
        can be a list of file paths (newick, nexus, hmatrix, 
        or a dictionary of trees.
    prefix : str
        Prefix of the file to save. If prefix is set to false, no file is
        saved. The default is False.
    taxa_replacement : str, optional
        Path of a table file containing two columns. The first column
        corresponds to the names of the terminals of the newick stored in
        infile, and the second column corresponds to their names the user wants
        to obtain at the end. All separators accepted. Example:
             AA Diceras
             AB Valletia
             AC Monopleura
        The default is False (no replacement).
    weighting : str
        Weighting scheme to use between: * FW (Fractional weighting from
                                               Rineau et al. 2021)
                                         * FWNL (Fractional weighting from
                                               Nelson and Ladiges),
                                         * UW (Uniform weighting from
                                               Nelson and Ladiges 1992),
                                         * MW (Minimal Weighting),
                                         * AW (Additive Weighting),
                                         * NW (No Weighting).
    parallel : str or int
        Option for choosing if the analysis is made using multiprocessing or
        not. This argument can be:
            * "no" if the user does not wan to use multiprocessing.
            * "auto" for automatic detection of the number of cpu available.
            * any integer corresponding to the number of cpu the user wants to
              allow to the analysis.
    dec_detail : bool
        If true, save a detailed output table in csv format named 
        prefix.triplet_table.csv of triplet weights per character in addition 
        to the .triplet file. 
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    verbose : bool
        Verbose mode if true.

    Returns
    -------
    triplet_dict : dict
        Dictionary of triplets (keys) and weights (values) without doublons.
    """

    if type(input_item) == list:
        character_dict = wrapper_character_extraction(input_item, 
                                                      taxa_replacement,
                                                      prefix,
                                                      verbose)
    
    else:
        character_dict = input_item

    # remove automatically repetitions if detected (user message printed)
    if rep_detector(character_dict):              
        character_dict = del_replications_forest(character_dict,
                                                     method=method,
                                                     prefix=prefix,
                                                     verbose=verbose)


    if not character_dict:
        print("ERROR: no informative character. Operation aborted.")
        sys.exit(1)

    # compute triplet dictionary for each tree in parallel
    print("Starting triplet decomposition and " + weighting + " computation")
    start = time.time()

    if parallel == "no":
        triplet_dict = standard_tripdec(character_dict, weighting, dec_detail, 
                                        prefix, verbose)
    else:
        triplet_dict = parallel_tripdec(character_dict, weighting,
                                        prefix, parallel)

    end = time.time()
    time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print("elapsed time (triplet decomposition and weighting): {}".format(
                                                                time_cptr))
    print("Generated " + str(len(triplet_dict))+" triplets")

    return triplet_dict
