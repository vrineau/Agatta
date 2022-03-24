# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@gmail.com

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
import warnings
import platform

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree


def search_pipeline(path_infile, software_path=False, software="paup",
                    prefix=False):
    """
    Function to run automatically a nexus file or a tnt file on PAUP, TNT, or
    WQFM respectivelly. The software to be used must be installed and
    accessible.

    Parameters
    ----------
    path_infile : string
        path of the nexus or tnt file.
    software_path : string
        path of the software chosen, paup, tnt, wqfm. Not necessari if paup is
        installed in a windows os.
    software : string, optional
        paup, tnt, wqfm. The default is "paup".

    Returns
    -------
    Nothing. Run the analysis with specified file and software

    """

    print("Running analysis on " + software + " software")
    start = time.time()

    if software_path:
        software_path = os.path.expanduser(software_path)
        if not os.path.isfile(software_path):
            print("ERROR: '" + software_path
                           + "' folder does not exist.\nOperation aborted.")
            sys.exit(1)

    if software == "tnt":
        print("WARNING: 3ia analysis using TNT is still under development.")
        if " " in path_infile:
            print("ERROR: TNT software doesn't allow paths with spaces.\n" +
                  "The path \"" + path_infile + "\" contains spaces.")
            sys.exit(1)

    ostype = platform.system()  # os detection

    if software == "wqfm":
        os.system(("java -jar \"" + software_path + "\" -i \"" +
                  path_infile + "\" -o \"" + prefix + ".tre\""))

    elif software == "tnt":
        os.system(software_path+" proc "+path_infile)

    elif software == "paup":
        if ostype == "Windows":
            begincmdline = "paup"
        else:
            begincmdline = software_path
        os.system(begincmdline + " -n \""+path_infile+"\" > /dev/null")

    print("Three-item analysis done")

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


def bandb(leaves, triplet_dict, base_tree=False,
          optimal_score=False, optimal_tree_list=[]):
    """
    Simple branch and bound analysis for three-item analysis. The function
    requires a list of leaves, a dictionary of triplets with their weights
    and performs a branch and bound analysis to found dichotomous trees that
    maximise the triplet score computed with the function tripletscore.

    Parameters
    ----------
    leaves : list
        List of leaf names.
    triplet_dict : dict
        Dictionary of triplets and weights coming from main_tripdec.

    Other parameters must be set by default for recursion.

    Returns
    -------
    int, float, or fraction
        Score of the optimal tree(s).
    list
        List of  optimal tree(s).

    """
    l = leaves.copy()

    if not optimal_score:
        optimal_tree = Tree()
        optimal_tree.populate(len(leaves), leaves)
        optimal_tree_list = [optimal_tree]
        optimal_score = tripletscore(triplet_dict, optimal_tree)
        base_tree = Tree()  # building two-leaves tree
        base_tree.add_child(name=l[0])
        base_tree.add_child(name=l[1])
        l.pop(0)
        l.pop(0)

    leaf = l[0]
    l.pop(0)

    # For each node, create an inserting node leading to it
    for node in base_tree.traverse():

        # Step 1 : rearrange
        if node.is_root():  # special case where the new leaf to branch is
                            # in sister-group of everything else
            newtree = Tree()
            insert = base_tree.copy()
            newtree.add_child(insert)
            newtree.add_child(name=leaf)

        else:  # other positions
            newtree = base_tree.copy()  # local tree

            if node.is_leaf():  # searching ancestral node
                newnode = newtree.get_leaves_by_name(node.name)[0]
            else:
                newnode = newtree.get_common_ancestor(
                    [l.name for l in node.get_leaves()])

            anc = newnode.get_ancestors()[0]  # for each node, take ancestor
            insertnode = anc.add_child(name="internal")  # branch node to it
            insert = newnode.copy()
            newnode.detach()  # del subtree
            insertnode.add_child(insert)  # connect the subtree on the leaf
            insertnode.add_child(name=leaf)  # add new leaf

        # Step 2 - compute score and recursion
        tripscore = tripletscore(triplet_dict, newtree)  # compute score
        newtree = [newtree]

        # Step 3 - check if score already too high and discard
        if not tripscore > optimal_score:  # if score not already too high

            # Step 4 - continue recursion to build a tree on all leaves
            if l:  # ifleaves are remaining
                tripscore, newtree = bandb(l,
                                     triplet_dict,
                                     newtree[0],
                                     optimal_score,
                                     optimal_tree_list)

            # Step 5 - recursion ends with an optimal tree
            if tripscore == optimal_score:
                add_tree = True
                for rftree in optimal_tree_list:
                    if rftree.robinson_foulds(newtree[0])[0] == 0:
                        add_tree = False
                if add_tree:
                    optimal_tree_list += newtree  # tree added to list

            elif tripscore < optimal_score:  # if recursion gives a result
                optimal_score = tripscore  # new optimal score
                optimal_tree_list = newtree  # new optimal tree

    return optimal_score, optimal_tree_list
