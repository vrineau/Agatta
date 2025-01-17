# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from random import choice
from .ini import taxa_extraction
import os
import sys
import time
import platform
from ete3 import Tree


def search_pipeline(path_infile, software_path=False, software="paup",
                    prefix=False):
    """
    Function to run automatically a nexus file, a tnt file, or a text file on
    PAUP, TNT, WQFM, or wTREE-QMC respectivelly. The software to be used must 
    be installed and accessible.

    Parameters
    ----------
    path_infile : string
        path of the input file.
    software_path : string
        path of the software chosen: paup, tnt, wqfm, wtree-qmc. Not necessary 
        if paup is installed in a windows os.
    software : string, optional
        paup, tnt, wqfm, wtree-qmc. The default is "paup".

    Returns
    -------
    Nothing. Run the analysis with specified file and software

    """

    print("Running analysis on " + software + " software")
    print("==================================================================")
    start = time.time()

    if software_path:
        software_path = os.path.expanduser(software_path)
        if not os.path.isfile(software_path):
            print("ERROR: '" + software_path
                           + "' folder does not exist.\nOperation aborted.")
            sys.exit(1)

    ostype = platform.system()  # os detection

    if ostype == "Windows":
        beginline = "paup"
        #endline = ">NUL"
    else:
        beginline = software_path
        endline = "> /dev/null"


    if software == "tnt":
        print("WARNING: 3ia analysis using TNT is still under development.\n"
              + "Searches using heuristics and fractional weights are not "
              + "optimized.")
        if " " in path_infile:
            print("ERROR: TNT software doesn't allow paths with spaces.\n" +
                  "The path \"" + path_infile + "\" contains spaces.")
            sys.exit(1)

        os.system(software_path + " proc "+path_infile + endline)

    if software == "wqfm":
        os.system("java -jar \"" + software_path + "\" -i \"" +
                  path_infile + "\" -o \"" + prefix + ".tre\"")

    elif software == "wtree-qmc":
        os.system("\"" + software_path + "\" -i \"" +
                  path_infile + "\" -o \"" + prefix + ".tre\"")

    elif software == "paup":

        os.system(beginline + " -n \"" + path_infile + "\"")

    print("==================================================================")
    print("Analysis done")

    end = time.time()
    time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print("elapsed time (3-item analysis): {}".format(time_cptr))


def rear_taxa(tree1):
    """
    Takes as argument an ete3 Tree and outputs the same tree with pruning
    and regrafting of a leaf.

    Parameters
    ----------
    tree1 : ete3.Tree
        Rooted tree.

    Returns
    -------
    tree2 : ete3.Tree
        Rearranged tree.
    leaf_node : ete3.Tree
        leaf node pruned and regrafted on the new tree

    """

    character_dict = {}
    character_dict[tree1.copy(method="newick")] = str(1)
    tree_copy_test = tree1.copy(method="newick")
    rf = 0  # similarity score between old and new tree

    while rf == 0:  # while trees are the same
        tree2 = tree_copy_test.copy(method="newick")
        rand_node = Tree()
        rand_taxa = ""

        # choice of pruned leaf and regraft location
        # check if pruned and regraft location are not the same
        while rand_node.name == rand_taxa:
            # location where to regraft the leaf
            rand_node = choice(
                [node for node in tree2.traverse(strategy="postorder")])
            # taxon randomly chosen to be  pruned and regrafted
            rand_taxa = choice(
                [taxa for taxa in taxa_extraction(character_dict)])

        regraft_node = rand_node.copy(method="newick")
        tree2.search_nodes(name=rand_taxa)[0].get_ancestors()[0].delete()
        tree2.search_nodes(name=rand_taxa)[0].delete()  # del taxon

        if regraft_node.search_nodes(name=rand_taxa):
            regraft_node.search_nodes(
                name=rand_taxa)[0].get_ancestors()[0].delete()
            regraft_node.search_nodes(name=rand_taxa)[0].delete()  # del taxon

        if len([leaf for leaf in rand_node.iter_leaves()]) == 1:
            leaf_name = rand_node.get_leaves()[0].name
            rand_node.get_leaves()[0].name = "new_ancestor"
            new_ancestor = rand_node.get_leaves()[0]
            new_ancestor.add_child(name=leaf_name)
            new_ancestor.add_child(name=rand_taxa)

        else:
            for node in rand_node.get_children():  # node deletion
                node.detach()
            new_ancestor = rand_node.add_child(name="new_ancestor")  # regraft
            new_ancestor.add_child(regraft_node)
            rand_node.add_child(name=rand_taxa)

        for node in tree2.traverse(strategy="levelorder"):
            if len(node.get_children()) == 1:
                node.children[0].delete()
        rf, max_rf, common_leaves, parts_t1, parts_t2, set1, set2 = \
            tree2.robinson_foulds(tree_copy_test)

    leaf_node = tree2.search_nodes(name=rand_taxa)[0]

    return tree2, leaf_node


def tripletscore(triplet_dict, stree):
    """
    Compute the triplet score of a dichotomic tree for a branch and bound
    analysis.
    The triplet score is the sum of weights of the triplets compatibles with
    the tree

    Parameters
    ----------
    triplet_dict : dict
        Dictionary of triplets generated by main_tripdec().
    stree : ete3.Tree
        Rooted tree.

    Returns
    -------
    score : int, float or fraction
        Score of the tree. The type depends of the type of the triplet weights.

    """

    score = 0

    for triplet, FW in triplet_dict.items():
        tripin = list(triplet.in_taxa)
        a = stree.get_leaves_by_name(tripin[0])
        b = stree.get_leaves_by_name(tripin[1])
        c = stree.get_leaves_by_name(list(triplet.out_taxa)[0])

        if a and b and c:  # if the three leaves are in the tree

            common = stree.get_common_ancestor([a[0], b[0]])
            cc = common.get_leaves_by_name(list(triplet.out_taxa)[0])

            if cc:  # if triplet not compatible
                score += FW*2
            else:  # if triplet compatible
                score += FW

    return score


def triplet_check(triplet, stree):
    """
    check if triplet exists in tree.
    Returns true or false.

    """

    tripin = list(triplet.in_taxa)
    a = stree.get_leaves_by_name(tripin[0])
    b = stree.get_leaves_by_name(tripin[1])
    c = stree.get_leaves_by_name(list(triplet.out_taxa)[0])

    if a and b and c:  # if the three leaves are in the tree

        common = stree.get_common_ancestor([a[0], b[0]])
        cc = common.get_leaves_by_name(list(triplet.out_taxa)[0])

        if cc:  # if triplet not compatible
            score = False
        else:  # if triplet compatible
            score = True

    return score

