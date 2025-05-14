# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from fractions import Fraction
from itertools import combinations
from .ini import taxa_extraction
from .ini import character_extraction
from .ini import wrapper_character_extraction
from .analysis import standard_tripdec
from .analysis import rep_detector
from .analysis import del_replications_forest
from .analysis import main_tripdec
import csv
from os import path, devnull, remove
import sys
from pypdf import PdfMerger
from fractions import Fraction
from ete3 import Tree

def constrict(treelist, prefix=False, silent=False):
    """
    Compute the strict consensus from a list of rooted trees.

    Parameters
    ----------
    treelist : list
        List of ete3.Tree objects.
    prefix : str or bool, optional
        Prefix of the file prefix.constrict to be saved with the strict
        consensus in newick format. The default is False (no file saved).

    Returns
    -------
    constrict : ete3.Tree
        Strict consensus.

    """

    if not silent:
        print("Strict consensus computation")

    if rep_detector(dict.fromkeys(treelist, 1)):
        print("ERROR: Repeated leaves have been detected.\n" +
                         "Operation aborted.")
        sys.exit(1)

    for t in treelist:
        for leaf in t.get_leaf_names():
            if not leaf in treelist[0].get_leaf_names():
                print("ERROR: All trees do not have the same " +
                         "length.\nOperation aborted.")
                sys.exit(1)


    constrict = treelist[0]  # base tree for the strict consensus computation
    del treelist[0]
    kill_node_list = []  # list of nodes to delete at the end

    if treelist:
        for node in constrict.traverse(strategy="preorder"):  # for each node
            if not node.is_leaf() and not node.is_root():  # check if exists in all
                for tree_compare in treelist:  # for each tree
                    node_find = False
                    for node_tree_compare in tree_compare.traverse(
                            strategy="preorder"):
                        set1 = set([leaf.name for leaf in node.get_leaves()])
                        set2 = set([leaf.name for leaf in
                                    node_tree_compare.get_leaves()])
                        if set1 == set2:
                            node_find = True
                            break
                    if node_find == False:
                        kill_node_list.append(node)
                        break

        for kill_node in kill_node_list:  # del nodes not present in all trees
            kill_node.delete()

    if prefix:
        with open(prefix+".constrict", "w") as constrictfile:
            constrictfile.write(constrict.write(format=9)+"\n")

    if not silent:
        print("Strict consensus computed")

    return constrict


def rcc(treelist, prefix=False, verbose=False):
    """
    Compute the reduced cladistic consensus (Wilkinson, 1994, 1995) from a 
    list of rooted trees.

          Wilkinson, M. (1994). Common cladistic information and its consensus
          representation: reduced Adams and reduced cladistic consensus trees
          and profiles. Systematic Biology, 43(3), 343-368.
          
          Wilkinson, M. (1995). More on Reduced Consensus Methods. Systematic 
          Biology, 44(3). 


    Parameters
    ----------
    treelist : list
        List of ete3.Tree objects.
    prefix : str or bool, optional
        Prefix of the file prefix.constrict to be saved with the strict
        consensus in newick format. The default is False (no file saved).
    verbose : bool, optional.
        Verbose mode if True. Default is False.

    Returns
    -------
    profile : list
        List containing the strict consensus and all trees from the profile.

    """

    def nts_to_trees(cardinal, in_list):
        """
        Compute a ete3 tree from a set of n-taxon statements

        """

        ntstree = Tree()
        ntstree.add_features(taxa_content = cardinal)

        for nts in sorted(in_list, key=lambda x: len(x), reverse=True):
            for node_iter in ntstree.traverse(strategy="postorder"):

                # set of taxa from all the direct descendant nodes
                taxa_set = set()
                for node_child1 in node_iter.get_children():
                    taxa_set = taxa_set | node_child1.taxa_content

                # if place to add child node, check all nts by size and add
                if len(node_iter.taxa_content - taxa_set) > 1 and (
                        (node_iter.taxa_content - taxa_set) >= nts or (
                            node_iter.taxa_content - taxa_set) == nts):
                    node_iter.add_child().add_features(taxa_content = nts)
                    break

        #label nodes: paralog, ortholog, apical
        for node in ntstree.traverse(strategy="preorder"):
            if [node1 for node1 in node.get_children() if node1.taxa_content]:
                taxa_set = set()
                for node_child2 in node.get_children():
                    taxa_set = taxa_set | node_child2.taxa_content

                for taxa in node.taxa_content - taxa_set:
                    node.add_child(name=taxa).add_features(
                        taxa_content = set())  # add leaf for each of the set

            elif node.taxa_content:  # if no child node
                for taxa in node.taxa_content:
                    node.add_child(name=taxa).add_features(
                        taxa_content = set()) # add leaf for each of the set

        return ntstree


    print("Reduced cladistic consensus (RCC) computation")

    if rep_detector(dict.fromkeys(treelist, 1)):
        print("ERROR: Repeated leaves have been detected.\n" +
                         "Operation aborted.")
        sys.exit(1)
    else:
        for t in treelist:
            for leaf in t.get_leaf_names():
                if not leaf in treelist[0].get_leaf_names():
                    print("ERROR: All trees do not have the same " +
                             "length.\nOperation aborted.")
                    sys.exit(1)

    # check size of trees
    cardinal_list = []
    for treel in treelist:
        leaf_cptr = 0
        for node in treel.traverse(strategy="preorder"):
            if node.is_leaf() == True:
                leaf_cptr += 1
        cardinal_list.append(leaf_cptr)

    cardinal_set = set(cardinal_list)

    if len(cardinal_set) > 1:
        print("ERROR: Trees must all have the same leaf set." +
                       "\nOperation aborted.")
        sys.exit(1)

    # build list of components per tree
    component_trees = []

    for treel in treelist:
        component_trees.append([])
        for node in treel.traverse(strategy="preorder"):
            if not node.is_leaf() and not node.is_root():
                component_trees[-1].append(
                    [set([leaf.name for leaf in node.get_leaves()]), set(
                        [leaf.name for leaf in treel.get_leaves()]) - set(
                            [leaf.name for leaf in node.get_leaves()])])

    # build intersection list
    nts_intersect = []

    for comp1 in component_trees[0]:
        for comp2 in component_trees[1]:
            if len(comp1[0] & comp2[0]) > 1 and len(comp1[1] & comp2[1]) > 0:
                nts_intersect.append([comp1[0] & comp2[0],comp1[1] & comp2[1]])

    del component_trees[1]  # del trees from component_trees after transfert
    del component_trees[0]  #idem

    #redundancy deletion
    del_nts_intersect = []  # list where to add nts to delete
    for nts1, nts2 in combinations(nts_intersect, 2):
        if (nts1[0] <= nts2[0] and (
                nts1[1] <= nts2[1] or nts1[1] == nts2[1])) and (
                nts1 not in del_nts_intersect):
            del_nts_intersect.append(nts1)

        elif (nts1[0] >= nts2[0] and (
                nts1[1] >= nts2[1] or nts1[1] == nts2[1])) and (
                nts2 not in del_nts_intersect):
            del_nts_intersect.append(nts2)

    for del_nts in del_nts_intersect:  # del redundant nts
        while del_nts in nts_intersect:
            nts_intersect.remove(del_nts)

    #intersection of each set of components per tree with nts_intersect
    while len(component_trees) > 0:
        nts_intersect_new = []
        for comp1 in component_trees[-1]:  # pile
            for comp2 in nts_intersect:
                if len(comp1[0] & comp2[0]) > 1 and len(
                        comp1[1] & comp2[1]) > 0:
                    nts_intersect_new.append(
                        [comp1[0] & comp2[0],comp1[1] & comp2[1]])

        nts_intersect = nts_intersect_new

        del component_trees[-1]

        #redundancy deletion
        del_nts_intersect = []
        doubles_nts_intersect = []
        for nts1, nts2 in combinations(nts_intersect, 2):
            if nts1 == nts2:
                doubles_nts_intersect.append(nts1)

            elif (nts1[0] <= nts2[0] and (
                    nts1[1] <= nts2[1] or nts1[1] == nts2[1])) and (
                    nts1 not in del_nts_intersect):
                del_nts_intersect.append(nts1)

            elif (nts1[0] >= nts2[0] and (
                    nts1[1] >= nts2[1] or nts1[1] == nts2[1])) and (
                        nts2 not in del_nts_intersect):
                del_nts_intersect.append(nts2)

        for del_nts in doubles_nts_intersect:
            while nts_intersect.count(del_nts) != 1:
                nts_intersect.remove(del_nts)

        for del_nts in del_nts_intersect:
            while del_nts in nts_intersect:
                nts_intersect.remove(del_nts)

        if len(nts_intersect) == 0:  # if nts_intersect empyt, star tree
            break

    # gather nts in lists with common cardinal
    nts_groups = {}
    for nts in nts_intersect:
        nts.append(frozenset(nts[0]|nts[1]))
        nts_groups[frozenset(nts[0]|nts[1])] = []

    for nts in nts_intersect:
        nts_groups[nts[2]].append(nts[0])

    for nts in nts_intersect:
        for cardinal, in_list in nts_groups.items():
            if nts[2] != cardinal and nts[2] >= cardinal and len(
                    nts[0] & set(cardinal)) > 1:
                nts_groups[cardinal].append(nts[0] & set(cardinal))
                
    # for each pair of trees: second step from Wilkinson 1995     
    for card1, card2 in combinations(nts_groups, 2):
        
        # cardinal computation
        newcard = frozenset(card1 & card2)
        ntslist = []

        # new nts list computation
        if len(newcard) > 3 and newcard != card1 and newcard != card2:
            for nts in nts_groups[card1] + nts_groups[card2]:
                
                newnts = set(nts & newcard)
                
                if newnts not in ntslist and len(newnts) > 1:
                    ntslist.append(newnts)
                    
        # check redundancy and validate
        if (len(ntslist) > len(nts_groups[card1])) and (
                len(ntslist) > len(nts_groups[card2])):
            nts_groups[newcard] = ntslist

    #merging all nts into trees
    profile = []

    for cardinal, in_list in nts_groups.items():  #for each tree
        tree_profile = nts_to_trees(cardinal, in_list)
        if not len(tree_profile) == leaf_cptr:
            profile.append(tree_profile)

    profile.append(constrict(treelist, prefix=False, silent=True))

    #sort trees by size
    def get_len(in_set):
        return len(in_set)

    profile = sorted(profile, key=get_len, reverse=True)  # sort by size

    if verbose:
        if not cardinal_list[0] in [len(profile_tree)
                                    for profile_tree in profile]:
            print("The strict consensus is not informative")

        #affichage des résultats
        if len(nts_intersect) == 0:
            print("RCC profile empty: no common phylogenetic information")
        else:
            for tree_profile in profile:
                print(tree_profile)

    if prefix:
        with open(prefix+".rcc", "w") as rccfile:
            first_line = True
            i = 1
            for profiletree in profile:
                if first_line:
                    rccfile.write("Strict consensus:    "
                                  + profiletree.write(format=9)+"\n")
                    first_line = False
                else:
                    rccfile.write("Profile tree " + str(i) + ":    "
                                  + profiletree.write(format=9) + "\n")
                    i += 1

    print("RCC of {} trees computed".format(str(len(profile))))

    return profile


def RI_path(cladopath, charpaths, taxarep1=False, taxarep2=False, 
       method="TMS", weighting="FW", prefix=False, verbose=False):
    """
    Compute the retention index for hierarchical characters (Kitching et al., 
    1998) in the three-item analysis framework from two files.
    The function write an output file in which is a global retention index
    of the analysis and a retention index for each hierarchical character tree.

          Kitching, I. J., Forey, P., Humphries, C., & Williams, D. (1998).
          Cladistics: the theory and practice of parsimony analysis.
          Oxford University Press.

    Parameters
    ----------
    cladopath : str
        Path to a newick file containing the cladogram. The tree must be the
        optimal cladogram obtained from the cladistic analysis of character
        trees stored in charpath.
    charpaths : list of str
        List of paths to files containing character trees in newick or 
        hmatrix format.
    taxarep1 : str, optional
        DESCRIPTION. The default is False.
    taxarep2 : str, optional
        DESCRIPTION. The default is False.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS" for removing polymorphism. 
        The default is "TMS".
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * UW (Uniform weighting from
                                       Nelson and Ladiges 1992),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        The default is "FW".
    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. The default is False (no file saved).

    Returns
    -------
    RI_char_dict : dict.
        Dictionary with character identifiers and global and per character RI.

    """

    print("Loading cladogram")

    cladopath = path.expanduser(cladopath)
    charpaths = [path.expanduser(charpath) for charpath in charpaths] 

    cladogram_dict = character_extraction(cladopath, taxarep1, verbose=False)

    print("Cladogram loaded")
    print("Loading character trees")

    character_dict = wrapper_character_extraction(charpaths, 
                                                  taxarep2,
                                                  prefix,
                                                  verbose=False)
        
    # Calls RI function
    RI_char_dict = RI(cladogram_dict, character_dict, taxarep1, taxarep2, 
           method, weighting, prefix, verbose)
    
    return RI_char_dict
        

def RI(cladogram_dict, character_dict, taxarep1=False, taxarep2=False, 
       method="TMS", weighting="FW", prefix=False, verbose=False):
    """
    Compute the retention index for hierarchical characters
    (Kitching et al., 1998) in the three-item analysis framework.
    The function write an output file in which is a global retention index
    of the analysis and a retention index for each hierarchical character tree.

          Kitching, I. J., Forey, P., Humphries, C., & Williams, D. (1998).
          Cladistics: the theory and practice of parsimony analysis.
          Oxford University Press.

    Parameters
    ----------
    cladogram_dict : dict
        Dictionary containing one newick tree (ete3 Tree objects) as key.
        The trees is generally the optimal tree or the strict consensus.
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
        The trees are generally the initial characters.
    taxarep1 : str, optional
        DESCRIPTION. The default is False.
    taxarep2 : str, optional
        DESCRIPTION. The default is False.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS" for removing polymorphism. 
        The default is "TMS".
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * UW (Uniform weighting from
                                       Nelson and Ladiges 1992),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        The default is "FW". Note that retention index by character state 
        in FW mode is equivalent to FWNL.
    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. The default is False (no file saved).

    Returns
    -------
    RI_char_dict : dict.
        Dictionary with character identifiers and global and per character RI.

    """

    # remove automatically repetitions if detected (user message printed)
    if rep_detector(character_dict):
        character_dict = del_replications_forest(character_dict,
                                                 method=method,
                                                 prefix=prefix,
                                                 verbose=verbose)

    print("Computing retention index")

    str_character_dict = {}
    for t, i in character_dict.items():
        str_character_dict[t] = str(i)
        
        for leaf in t.get_leaf_names():
            if not leaf in list(cladogram_dict.keys())[0].get_leaf_names():
                print("ERROR: Character tree " + str(i)
                               + ": leaf set must be equal"
                               + " or a subset of the cladogram\n"
                               + "leaf set. Operation aborted.")
                sys.exit(1)

    c_triplet_dict = standard_tripdec(cladogram_dict,
                                      weighting,
                                      dec_detail=False,
                                      prefix=False,
                                      verbose=False)
    RI_tot = [0, 0]
    RI_char_dict = {}
    RI_char_dict_num = {}
    RI_char_dict_denom = {}

    # computation of the global RI in the same time as RI per character
    for chartree, keys in str_character_dict.items():
        RI_char_dict[keys] = 0
        RI_char_dict_num[keys] = 0
        RI_char_dict_denom[keys] = 0
        triplet_dict = standard_tripdec({chartree: keys},
                                        weighting,
                                        dec_detail=False,
                                        prefix=False,
                                        verbose=False)

        for triplet1, FW in triplet_dict.items():  # for each triplet in a char
            if triplet1 in c_triplet_dict:  # if triplet is in cladogram, add

                #RI numerator
                RI_char_dict_num[keys] += FW  # per character
                RI_tot[0] += FW  # total RI

            #RI denominator
            RI_char_dict_denom[keys] += FW
            RI_tot[1] += FW


    # computation of RI per character state
    # This functionality is only able when using hmatrix, which gives specific 
    # codes for each character state.
    
    for chartree, keys in str_character_dict.items():
        
        for node in chartree.traverse(strategy="preorder"):
            if node.is_leaf() == False and node.is_root() == False:
                try: 
                    charstate_count = node.charstate_name

                    # compute triplets for this single state
                    charstate = chartree.copy('cpickle')
                    for n in charstate.traverse(strategy="preorder"):
                        if n.is_leaf() == False and n.is_root() == False:
                            try:
                                if not n.charstate_name == charstate_count:
                                    n.delete()
                            except AttributeError:
                                n.delete()
                                
                    keystate = str(keys.split('.')[0]) + "_" + str(
                                                            charstate_count)
                    RI_char_dict[keystate] = 0
                    RI_char_dict_num[keystate] = 0
                    RI_char_dict_denom[keystate] = 0
                    
                    state_triplet_dict = standard_tripdec({charstate: keys},
                                                    weighting,
                                                    dec_detail=False,
                                                    prefix=False,
                                                    verbose=False)

                    for triplet1, FW in state_triplet_dict.items():  
                        if triplet1 in c_triplet_dict:
            
                            #RI numerator
                            RI_char_dict_num[keystate] += FW  # per state
            
                        #RI denominator
                        RI_char_dict_denom[keystate] += FW
            
                except AttributeError:
                    pass
    
    # merge all in RI_char_dict fractions dict and remove zero fractions
    for keys, values in RI_char_dict_num.items():
        if RI_char_dict_denom[keys] == 0:
            RI_char_dict[keys] = [0, 1]
        else:
            RI_char_dict[keys] = [RI_char_dict_num[keys],
                                    RI_char_dict_denom[keys]]

    # computation of ri when free-paralogy subtrees (poly)
    for keys, values in str_character_dict.items(): # build RI for entire chars
        if '.' in str(values):
            if not values.split('.')[0] in RI_char_dict:
                RI_char_dict[values.split('.')[0]] = [0, 0]

            RI_char_dict[values.split('.')[0]][0] += RI_char_dict[values][0]
            RI_char_dict[values.split('.')[0]][1] += RI_char_dict[values][1]

    if RI_tot[1] == 0:
        RI_char_dict["Total"] = [0, 1]  # add global RI
    else:
        RI_char_dict["Total"] = [RI_tot[0], RI_tot[1]]

    print("Retention index computed")

    # computation of ri per state (results in FW are not equivalent between 
    # states and characters because of the correction (Rineau et al. 2021).
    #
    # Rineau, V., Zaragüeta, R., & Bardin, J. (2021). Information content 
    # of trees: three-taxon statements, inference rules and dependency. 
    # Biological Journal of the Linnean Society, 133(4), 1152-1170.) 

    # Write file
    if prefix:
        with open(prefix+".ri", "w") as RI_file:
            
            def formatNumber(num): # remove zero after float
                if num % 1 == 0:
                    return int(num)
                else:
                    return num

            def writeri(keys, values):
                RI_string0 = "[" + str(keys) + "]"
                RI_string1 = str(formatNumber(round((float(Fraction(values[0],
                                                       values[1]))*100), 2)))
                RI_string2 = str(formatNumber(round((float(values[0])), 2)))
                RI_string3 = str(formatNumber(round((float(values[1])), 2)))
                
                RI_string = [RI_string0, RI_string1, RI_string2, RI_string3]
                
                print('\t'.join(RI_string).expandtabs(10))
                RI_file.write((' '.join(RI_string) + "\n").expandtabs(10))
            
            RI_string = ['Chars', 'RI', 'Retained', 'Total']
            
            print('\n' + '\t'.join(RI_string).expandtabs(10))
            RI_file.write((' '.join(RI_string) + "\n").expandtabs(10))
            
            for keys, values in RI_char_dict.items():  # IR characters
                if not '.' in keys and not '_' in keys:
                    writeri(keys, values)
                    
            RI_string = ['\nStates', 'RI', 'Retained', 'Total']
            
            print('\t'.join(RI_string).expandtabs(10))
            RI_file.write((' '.join(RI_string) + "\n").expandtabs(10))
                    
            statescodespresence = False
            for keys, values in RI_char_dict.items():  # IR states
                if not '.' in keys and '_' in keys:
                    writeri(keys, values)
                    statescodespresence = True
                    
            if not statescodespresence:
                info_string = ("NA (hmatrix mandatory to set states codes)")
                print(info_string)
                RI_file.write(info_string + "\n")

            RI_string = ['\nSubtrees', 'RI', 'Retained', 'Total']
            
            print('\t'.join(RI_string).expandtabs(10))
            RI_file.write((' '.join(RI_string) + "\n").expandtabs(10))

            polypresence = False
            for keys, values in RI_char_dict.items():  # IR subtrees
                if '.' in keys and not '_' in keys:
                    writeri(keys, values)
                    polypresence = True

            if not polypresence:
                info_string = "NA\n"
                print(info_string)
                RI_file.write(info_string)
                
    # Codes for RI_char_dict
    #
    # 1 -   character
    # 1.1 - character subtree coming from FPLSA
    # 1_1 - character 1, state 1
    #
    
    RI_char_dict_return = {}
    for keys, values in RI_char_dict.items():
        RI_char_dict_return[keys] = Fraction(values[0], values[1])
    
    return RI_char_dict_return


def triplet_distance(t1, t2, prefix=False, method="TMS", weighting="FW", 
                     silent=False, verbose=False):
    """
    Compute several metrics between two hierarchical trees including the 
    inter-tree retention index (ITRI; Grand et al. 2014, Rineau et al. 2015, 
    which is an asymmetric distance between trees) and the triplet distance.
    Can be used also in the specific case of comparing a reference tree to 
    another tree.

    Grand, A., Zaragüeta-Bagils, R., Velez, L. M., & Ung, V. (2014).
    A cladistic re-analysis of the Gadiformes (Teleostei, Paracanthopterygii)
    using three-item analysis. Zootaxa, 3889(4), 525-552.
    Rineau, V., Grand, A., Zaragüeta i Bagils, R., & Laurin, M. (2015). 
    Experimental systematics : Sensitivity of cladistic methods to 
    polarization and character ordering schemes. 
    Contributions to Zoology, 84(2), 129‑148.


    Parameters
    ----------
    t1 : dict
        Dictionary containing one newick tree (ete3 Tree object) as key.
        If one want comparing a tree to a reference tree, t1 is the 
        reference tree.
    t2 : dict
        Dictionary containing one newick tree (ete3 Tree object) as key.
        If one want comparing a tree to a reference tree, t2 is the 
        tree one wants to compare to the reference.
    prefix : str
        Prefix of the text file to be saved containing the resulting ITRI.
        The complete path can be used. The default is False (no file saved).
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS" for removing polymorphism. 
        The default is "TMS".
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * UW (Uniform weighting from
                                       Nelson and Ladiges 1992),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        The default is "FW".
    Returns
    -------
    float
        precision. ITRI(t1,t2) (Precision: (t1∩t2)/t2)) 
        Computing the ITRI{t1,t2}, i.e. proportion of the FW of the  
        3is of t1 that are also present in t2. If t1 is a reference tree, 
        ITRI12 is the Precision: proportion of triplets from the reconstructed 
        tree that are true (equivalent to power in Rineau et al. 2015).
    float
        recall. ITRI(t2,t1) (Recall: (t1∩t2)/t1))
        Computing the ITRI{t1,t2}, i.e. proportion of the FW of the 3is 
        of t2 that are also present in t1. If t1 is a reference tree, ITRI12 is
        the recall.: proportion of triplets from the true tree that have been
        retrieved (inverse of artefact in Rineau et al. 2015).
    float
        tripdistance. Harmonic mean of ITRI(t2,t1) and ITRI(t1,t2).
        Harmonic mean between the two asymmetric ITRI values
        computed. Gives a triplet distance between two trees.
        If t1 is a reference tree, it is the F-score (harmonic mean of 
        precision and recall).

    """
    
    def triplet_incongruence_detector(triplet, tripdic):
        "Test if a triplet is incompatible with the triplet set (True), or"
        "absent but compatible (False). The triplet must not be in the set"
        
        for t in tripdic:
            t1set = t.out_taxa.union(t.in_taxa)
            t2set = triplet.out_taxa.union(triplet.in_taxa)
            if t1set == t2set:
                return True
            
        return False

    # remove automatically repetitions if detected (user message printed)
    if rep_detector({t1: 0}):
        t1 = del_replications_forest({t1: 0},
                                                 method=method,
                                                 prefix=prefix,
                                                 verbose=verbose).items()[0]

    if rep_detector({t2: 0}):
        t2 = del_replications_forest({t2: 0},
                                                 method=method,
                                                 prefix=prefix,
                                                 verbose=verbose).items()[0]

    tt_tripdic = standard_tripdec({t1: 0},
                                      weighting,
                                      dec_detail=False,
                                      prefix=False,
                                      verbose=False)

    rt_tripdic = standard_tripdec({t2: 0},
                                      weighting,
                                      dec_detail=False,
                                      prefix=False,
                                      verbose=False)
    

            

    RI_reconstructed_tree = Fraction(0, 1)  # total reconstruted tree
    RI_true_tree = Fraction(0, 1)  # total true tree
    RI_intersect_reconstructed_tree = Fraction(0, 1)
    RI_intersect_true_tree = Fraction(0, 1)
    false_positives = Fraction(0, 1)
    false_negatives = Fraction(0, 1)
    incongruent_triplets_t1 = Fraction(0, 1)
    incongruent_triplets_t2 = Fraction(0, 1)
    neutral_triplets_t1 = Fraction(0, 1)
    neutral_triplets_t2 = Fraction(0, 1)
    
    for trip, FW in rt_tripdic.items():
        if trip in tt_tripdic:
            RI_intersect_reconstructed_tree += FW
        else:
            if triplet_incongruence_detector(trip, tt_tripdic):
                incongruent_triplets_t1 += FW
            else:
                neutral_triplets_t1 += FW
            false_positives += FW
        RI_reconstructed_tree += FW

    for trip, FW in tt_tripdic.items():
        if trip in rt_tripdic:
            RI_intersect_true_tree += FW
        else:
            if triplet_incongruence_detector(trip, rt_tripdic):
                incongruent_triplets_t2 += FW
            else:
                neutral_triplets_t2 += FW
            false_negatives += FW
        RI_true_tree += FW

    recall = float(Fraction(RI_intersect_true_tree, RI_true_tree))  # power

    try:
        precision = float(Fraction(RI_intersect_reconstructed_tree,
                               RI_reconstructed_tree))  # 1/artefact
    except ZeroDivisionError:
        precision  = 0

    # F1-score
    try:
        F1score = (2 * precision * recall) / (precision + recall)
    except ZeroDivisionError:
        F1score  = 0
        
    #Triplet distance sensu Dobson 1975 Dobson, A.J. (1975). Comparing the 
    # shapes of trees. In: Street, A.P., Wallis, W.D. (eds) Combinatorial 
    # Mathematics III. Lecture Notes in Mathematics, vol 452.
    tripdistance = abs(false_positives + false_negatives) / 2
    
    # Triplet congruence / similarity
    # Amount of common triplets between t1 and t2. It may be a 
    # different amount in t1 and t2 because of the fractionnal weighting
    tripcongruence = abs(RI_intersect_true_tree + 
                         RI_intersect_reconstructed_tree) / 2
    
    # Triplet congruence / similarity. Not the same as triplet distance, as
    # the triplet incongruence takes only into account the incompatible 
    # triplets, while triplet distance takes also into account triplet that are
    # different but still compatible (or more exactly, neutral)
    tripincongruence = incongruent_triplets_t1 + incongruent_triplets_t2

    tripneutral = neutral_triplets_t1 + neutral_triplets_t2
        
    # results list to display /save
    itristr = []

    itristr.append("Number of triplets in t1 (relevant elements) : "
                   + str(round(float(RI_true_tree), 3)))
    
    itristr.append("Number of triplets in t2 (retreived elements) : " 
                   + str(round(float(RI_reconstructed_tree), 3)))
    
    itristr.append("Number of triplets both in t1 and t2 " 
                   + "(true positives) : " 
                   + str(round(float(RI_intersect_true_tree), 3)) 
                   + " (t1 also in t2) "
                   + str(round(float(RI_intersect_reconstructed_tree), 3)) 
                   + " (t2 also in t1)")
    
    itristr.append("Number of triplets in t2 but not in t1 " 
                   + "(false positives) : " 
                   + str(round(float(false_positives), 3)))
    
    itristr.append("Number of triplets in t1 but not in t2 " 
                   + "(false negatives) : " 
                   + str(round(float(false_negatives), 3)))
    
    itristr.append("Triplet distance ((t2Δt1)/2)) : " 
                   + str(round(float(tripdistance), 3)))
    
    itristr.append("Triplet congruence (amount of triplets present in both" 
                   + " trees) (t2 in t1 + t1 in t2)/2)) : " 
                   + str(round(float(tripcongruence), 3)))

    itristr.append("Triplet incongruence (amount of triplets incompatibles " +
                   "between t1 and t2) : " 
                   + str(round(float(tripincongruence), 3)))

    itristr.append("Triplet neutral (amount of neutral triplets, i.e. " 
                   + "triplets compatibles which are not the same between " 
                   + "t1 and t2) : " 
                   + str(round(float(tripneutral), 3)))

    itristr.append("ITRI(t1,t2) (Precision: (t2 in t1)/t2)) : " 
                   + str(round(precision, 3)))
    
    itristr.append("ITRI(t2,t1) (Recall: (t1 in t2)/t1)) : " 
                   + str(round(recall, 3)))
    
    itristr.append("Triplet F1-score: " 
                   + "(2 * Precision * Recall) / " 
                   + "(Precision + Recall)) : " 
                   + str(round(F1score, 3))) 

    # save file
    if prefix:
        with open(prefix+".itri", "w") as itrifile:
            
            for line in itristr:
                itrifile.write(line + "\n")

    # display
    if not silent:
        for line in itristr:
            print(line)
            
    return (float(RI_true_tree), 
           float(RI_reconstructed_tree),
           float(RI_intersect_true_tree),
           float(RI_intersect_reconstructed_tree),
           float(false_positives), 
           float(false_negatives),
           float(tripdistance),  #Dobson 1975
           float(tripcongruence),
           float(tripincongruence),
           float(tripneutral),
           precision,  # ITRI(T1,T2)
           recall,     # ITRI(T2,T1)
           F1score)              
           

def cladogram_label(cladogram, clade_number_option, clade_type_option, 
                    pdf_files=False):
    """
    Label tree nodes for NRI and character state test.
    """
    
    #sometimes issues with pyqt installation
    if pdf_files:
        try:
            from ete3 import NodeStyle, TreeStyle, faces
            from ete3 import TextFace, COLOR_SCHEMES, CircleFace
            
        except ImportError:  # issue with ete3 imports
            print("A manual instal of PyQt5 is requested to use the --pdf "
                  + "functionality\nPlease install using 'pip install" +
                  " PyQt5' and rerun the command line")
            sys.exit(1)

    cladogram.ladderize()

    if clade_number_option == "yes":  # label cladogram nodes
        node_style_num_count = 1
        for node in cladogram.traverse(strategy="preorder"):
            if node.is_leaf() == False  and node.is_root() == False:
                if pdf_files:
                    style_num = TextFace(node_style_num_count, fsize=2)
                    node.add_face(style_num, column=1,
                                  position = "branch-bottom")
                node.add_feature("clade_label", node_style_num_count)
                node.name = node_style_num_count
                node_style_num_count += 1
            elif node.is_root() == True:
                node.add_feature("clade_label", "root")
            else:
                node.add_feature("clade_label", node.name)

    if clade_type_option == "yes":  # label cladogram node types
        for node in cladogram.traverse(strategy="postorder"):
            if node.is_leaf():
                node.add_feature("clade_type", "leaf")
            else:
                i = len([child_node for child_node in node.get_children()
                         if not child_node.is_leaf()])
                if i >= 2:
                    node.add_feature("clade_type", "paralog")
                elif i == 1:
                    node.add_feature("clade_type", "ortholog")
                elif i == 0:
                    node.add_feature("clade_type", "apical_node")

    syn_dict = {}
    for node_number in cladogram.traverse(strategy="postorder"):
        syn_dict[str(node_number.clade_label)] = {
            "accepted": [],"rejected": []}

    return cladogram, syn_dict


def NRI(cladopath, charpaths, taxarep1=False, taxarep2=False, prefix="rnri", 
        rnri_codes=False, weighting='FW', polymethod='TMS', totaltree=True,
        rescaled=True, pdf_files=False):
    """
    Compute the Nodal Retention Index and return a pdf with percentages of NRI
    for each character (or for each set of characters if rnri_codes=True) 
    plus a csv file and with the total amount of triplet weights by node.
    
    Nodal retention index: given a set of characters and a cladogram, gives the 
    amount of triplet from each character in agreement with each node. Allows 
    to see easily the relative support for each clade.

    Parameters
    ----------
    cladopath : str
        path of input file containing the cladogram (or consensus).
    charpaths : list of str
        list of path of input file containing characters in hmatrix or newick 
        format.
    taxarep1 : str, optional
        DESCRIPTION. The default is False.
    taxarep2 : str, optional
        DESCRIPTION. The default is False.
    prefix : str
        Prefix of the files prefix.nri (csv file + a newick tree with node 
        names), prefix.pdf (percentages) and prefix.total.pdf (total triplet
        weights) to be saved. 
        The default is 'rnri'.
    rnri_codes : bool
        By default, the prefix.pdf file will display a tree with a piechart for 
        each node representing the relative proportion of support from each 
        character. The user can 
    weighting : str, optional
        Weighting scheme to use between FW, FWNL, MW, AW, NW. 
        The default is 'FW'.
    polymethod : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS". The default is "TMS".
    totaltree : bool, optional
        If set to True, computes the prefix.total.pdf showing the support of 
        each node given the relative amounts of character triplets in 
        agreement. The default is True.
    rescaled : bool, optional
        If True, for prefix.total.pdf the absolute values will be be rescaled 
        to avoid the symetry effect inherent to triplet decomposition. 
        Because of this effect, a balanced character has a heavier weight
        than an inbalanced character.
        The default is True.

    Returns
    -------
    pie_percentages1 : dict
        dictionary with node names as keys and lists of percentages as values.
        The list correspond all percentages displayed in the piechart tree.

    """

    #sometimes issues with pyqt installation
    if pdf_files:
        try:
            from ete3 import NodeStyle, TreeStyle, faces
            from ete3 import TextFace, COLOR_SCHEMES, CircleFace
            
        except ImportError:  # issue with ete3 imports
            print("A manual instal of PyQt5 is requested to use the --pdf "
                  + "functionality\nPlease install using 'pip install" +
                  " PyQt5' and rerun the command line")
            sys.exit(1)
    
    # reading matrix cladogram csv
        
    character_dict = wrapper_character_extraction(charpaths, taxarep1, prefix,
                                                  verbose=False)
        
    cladogram_dict = character_extraction(cladopath, taxarep2, verbose=False)
    
    if rnri_codes:
        
        rnri_codes = path.expanduser(rnri_codes)
        
        with open(rnri_codes, "r") as rnri_file:
            try:
                dialect = csv.Sniffer().sniff(rnri_file.readline(),
                                              delimiters=[";","\t"," ","|"])
            except:
                print("ERROR: The table file '" + rnri_codes +
                                  "' is not adequately formated." +
                                  "\n Operation aborted.")
            else:
                no_exception = True
            if not 'no_exception' in locals():
                sys.exit(1)
        
        with open(rnri_codes, "r") as rnri_file:
        
            char_list_raw = csv.reader(rnri_file, delimiter=dialect.delimiter)
            rnri_codes = list(char_list_raw)
        
            if not len({len(l) for l in rnri_codes}) == 1:
                 print("ERROR: The table file '" + rnri_codes +
                                  "' is broken.\n Operation aborted.")
                 sys.exit(1)
    
    # default function use : each character is represented in the piechart
    else:
        rnri_codes = []
        for charname in character_dict.values():
            rnri_codes.append(['', 'character ' + str(charname)])
            
    # add names to internal nodes
    cladogram = list(cladogram_dict.keys())[0]
    cladogram.ladderize()
    pie_percentages = dict() # dict int node name / percentage list for piech.
    Node_weight = dict() # dict internal node name / percentage list for piech.
    Node_weight_rescaled = dict() # dict with rescaling denominator
    
    categories_rnri = []
    n = [categories_rnri.append(x) for x in [x[1] for x in rnri_codes] 
     if x not in categories_rnri]
    
    node_style_num_count = 1
    
    card = len(cladogram.get_leaf_names())
    rescaledtot = 0
    for node in cladogram.traverse(strategy="preorder"):
        if node.is_leaf() == False  and node.is_root() == False:
            cardn = len(node.get_leaf_names())
            rescaledtot += (cardn - 1)*(card - cardn)
        
    for node in cladogram.traverse(strategy="preorder"):
        if node.is_leaf() == False  and node.is_root() == False:
            node_style_num_countstr = str(node_style_num_count)
            
            if pdf_files:
                style_num = TextFace(node_style_num_countstr, fsize=2)            
                node.add_face(style_num, column=1,
                              position = "branch-bottom")
                
            node.add_feature("clade_label", node_style_num_countstr)
            node.name = node_style_num_countstr
            pie_percentages[node_style_num_countstr] = [
                Fraction(0,1)] * len(set([x[1] for x in rnri_codes]))
            Node_weight[node_style_num_countstr] = Fraction(0,1)
            node_style_num_count += 1
            
            # rescaling factor computation
            cardn = len(node.get_leaf_names())
            Node_weight_rescaled[
                node.name] = ((cardn - 1)*(card - cardn)) / rescaledtot
            
        else:
            node.add_feature("clade_label", False)
    
    # RNRI computation : decompose each character in triplets with weights
    # for each triplet : traverse the tree and check if agreement tree/triplet
    # if yes add value to pie_percentages
    
    print("Starting detailed triplet decomposition and FW computation for NRI")
    
    for character, carnb in character_dict.items():  # for each character
    
        print("Decomposition: character " + str(carnb))
        
        #blockPrint()
        triplist = main_tripdec({character : 1}, prefix=False, 
                                taxa_replacement=False, weighting=weighting, 
                                parallel='no', dec_detail=False, 
                                method=polymethod, verbose=False)
        #enablePrint()
            
        for triplet_rnri, FW in triplist.items():  # check if triplet is true
            
            in1name, in2name = list(triplet_rnri.in_taxa)
            
            in1 = cladogram.get_leaves_by_name(str(in1name))[0]
            in2 = cladogram.get_leaves_by_name(str(in2name))[0] 
            out = cladogram.get_leaves_by_name(list(
                triplet_rnri.out_taxa)[0])[0]
            
            in_node = set(Tree.get_leaf_names(in1.get_common_ancestor(in2)))
            nodeout1 = set(Tree.get_leaf_names(out.get_common_ancestor(in1)))
            nodeout2 = set(Tree.get_leaf_names(out.get_common_ancestor(in2)))
            
            # check if triplet is valid
            if (in_node.issubset(nodeout1) and in_node.issubset(nodeout2) and 
                in_node != nodeout1 and in_node != nodeout2):
                
                # add value in result list
                pie_percentages[in1.get_common_ancestor(in2).name][
                    categories_rnri.index(rnri_codes[int(carnb)-1][1])] += FW
                Node_weight[in1.get_common_ancestor(in2).name] += FW
    
    print("Triplet decomposition and FW computation completed")
    
    # round and compute percentages
    
    pie_percentages1 = dict()
    Node_weight2 = dict()
    
    for node_name, values_list in pie_percentages.items():
        
        total = sum(pie_percentages[node_name])
        pie_percentages1[node_name] = [
            round((x / total) * 100) for x in values_list]
        pie_percentages1[node_name].pop()
        totalsub = sum(pie_percentages1[node_name])
        pie_percentages1[node_name].append(100 - totalsub)
    
    for node_name, frac in Node_weight.items():
        Node_weight2[node_name] = round((float(frac)), 5)
    
    # display pie-percentages
    
    # Basic tree style
    if pdf_files:
    
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.margin_top = 10
        ts.margin_right = 10
        ts.margin_left = 10
        ts.margin_bottom = 10
        ts.show_scale = False
        colors=["#a6cee3",
                "#b2df8a",
                "#fb9a99",
                "#fdbf6f",
                "#cab2d6",
                "#ffff99",
                "#1f78b4",
                "#33a02c",
                "#e31a1c",
                "#ff7f00",
                "#6a3d9a",
                "#b15928",
                "#001b17",
                "#9dff47",
                "#794ffd",
                "#009243",
                "#01168f",
                "#a6ffb7",
                "#c2008e",
                "#044a00",
                "#af7eff",
                "#ce6500",
                "#0091b6",
                "#c71b00",
                "#000c2f",
                "#f7ffde",
                "#3e0012",
                "#ffcdc4",
                "#b40021",
                "#4b4800",
                "#a60052",
                "#986400"]
        
        #remove display nodes
        for n in cladogram.traverse():
           nstyle = NodeStyle()
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 0
           n.set_style(nstyle)
           
        #Associate the PieChartFace only with internal nodes
        def pie_nodes(node):
            if not node.is_leaf() and not node.is_root():
                
                F= faces.PieChartFace(pie_percentages1[node.name],
                                      colors=colors, 
                                      width=50, height=50)
                F.border.width = None
                F.opacity = 0.6
                faces.add_face_to_node(F,node, 0, position="branch-right")
        
        ts.layout_fn = pie_nodes
        
        i = 0
        for cat in categories_rnri:
            if i > len(colors):
                print("WARNING: too many categories, "
                      +"piechart will be unreadable.")
            else:
                ts.legend.add_face(CircleFace(10, colors[i]), column=0)
                ts.legend.add_face(TextFace(cat, fsize=7), column=1)    
            
                i += 1 
            
        # if yes, save a pdf with node size reflecting overall triplet support
        if totaltree:
            
            ts2 = TreeStyle()
            ts2.show_leaf_name = True
            ts2.margin_top = 10
            ts2.margin_right = 10
            ts2.margin_left = 10
            ts2.margin_bottom = 10
            ts2.show_scale = False
            
            bubble_cladogram = cladogram.copy()
            
            def bubble_layout(node):
                if not node.is_leaf() and not node.is_root():
                    # Creates a sphere face whose size is proportion. to weight
                    C = CircleFace(radius=radius[node.name], 
                                   color="#6d87ed", style="sphere")
                    C.opacity = 0.6
                    faces.add_face_to_node(C, node, 0, position="float")
                    
            #réechelonner
            radius = dict()
            if rescaled:
                radiusr = dict()
                for nodename, val in Node_weight2.items():
                    radiusr[nodename] = (Node_weight2[nodename] / 
                                         Node_weight_rescaled[nodename])
                for nodename, val in radiusr.items():
                    radius[nodename] = (radiusr[nodename] / 
                                        max(radiusr.values()))*20
                    
            else:
                for nodename, val in Node_weight2.items():
                    radius[nodename] = (Node_weight2[nodename] / 
                                        max(Node_weight2.values()))*20
                    
            ts2.layout_fn = bubble_layout
            
            bubble_cladogram.render(prefix + ".total.pdf", tree_style=ts2)
        
        # save results
        cladogram.render(prefix + ".pdf", tree_style=ts)
    
    with open(prefix+".nri", "w") as results_file:
        
        results_file.write(' ')
        for e in categories_rnri:
            results_file.write(',' + e)
        results_file.write(",total,total_rescaled\n")  
        
        for carname, percentage_list in pie_percentages1.items():
            results_file.write('Node_' + carname)
            for p in percentage_list: 
                results_file.write(',' + str(p))
    
            radius = Node_weight2[carname]        
            radius_rescaled = Node_weight2[carname] / Node_weight_rescaled[
                carname]
                
            results_file.write("," + str(radius) + "," + 
                               str(radius_rescaled) + "\n")  
    
            
        results_file.write("\n" + cladogram.write(format=8) + "\n")
            
    print("Nodal Retention Index computation ended successfully")
        
    return pie_percentages1
    

def character_states_test(cladogram_dict, character_dict,
                          prefix, pdf_files=False):
    """
    Run the character states testing procedure for hierarchical characters
    defined by Cao (2008) and corrected by Rineau (2017). This procedure tests
    at each symmetric node (paralog) of the cladogram the presence of
    instances of the derived state.

          Cao, N. (2008). Analyse à trois éléments et anatomie du bois des
          Fagales Engl (Doctoral dissertation, Paris, Muséum national
          d'histoire naturelle).

          Rineau, V. (2017). Un nouveau regard cladistique sur l'anatomie
          comparée, la phylogénie, la systématique et la paléoécologie des
          rudistes (Bivalvia, Hippuritida) (Doctoral dissertation,
          Université Pierre et Marie Curie).

    Parameters
    ----------
    cladogram_dict : dict
        Dictionary containing one newick tree (ete3 Tree objects) as key.
        The tree must be the optimal cladogram obtained from the cladistic
        analysis of character trees stored in the character_dict argument.
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
        The trees must be the initial characters.
    prefix : str
        Prefix of the files to be computed.
    pdf_files : bool, optional
        If true, save a pdf file with one page for each character state for
        visualizing the results of the character state procedure.
        The default is False (no pdf files computed).

    Returns
    -------
    dict
        Two dictionaries with the results of the procedure sorted in two ways.

    """

    def find_synapomorphy(leaf_syn):
        """
        For one leaf, recursion over the tree in preorder and find the
        corresponding synapomorphy.

        """

        if leaf_syn.is_root():
            ancestor_syn = leaf_syn
        elif leaf_syn.get_ancestors()[0].state == "y":
            ancestor_syn = find_synapomorphy(leaf_syn.get_ancestors()[0])
        else:  # when acestor is find, captured here
            ancestor_syn = leaf_syn
        return ancestor_syn


    def synapomorphy_test_placement(cladogram, taxa_in_state,
                                    results_test_dict, char_names,
                                    char_states_names):
        """
            Test and placement of a character state.

        """
        results_test_dict[char_names][char_states_names] = []
        synapomorphies_set = set()  # set of nodes with derived state
        synapo_test = str("accepted")  #accepted(default)/rejected
        potential_synapo = cladogram.get_common_ancestor(taxa_in_state)

        # assignation of the state to a node
        # if no incongruence or no symmetric node implied
        if potential_synapo.get_leaf_names() == taxa_in_state or not (
                potential_synapo.search_nodes(clade_type="paralog")):
            synapomorphies_set.add(potential_synapo)

        # if ambiguous interpretation, check leaves in the symmetric nodes
        else:
            for node_cladogram in cladogram.traverse(strategy="postorder"):
                node_cladogram.add_feature("state", "n")  # no

            for leaf_node1 in cladogram.iter_leaves():
                if leaf_node1.name in taxa_in_state:
                    leaf_node1.state = "y"  # yes

            for node_cladogram in cladogram.traverse(strategy="postorder"):

                if node_cladogram.clade_type == "apical_node":
                    if len([apical_node_leaf for apical_node_leaf in
                            node_cladogram.children if
                            apical_node_leaf.state == "y"]) > 0:
                        node_cladogram.state = "y"

                elif node_cladogram.clade_type == "paralog":
                    test_node_cptr = 0
                    for test_node in (gen_paralog_node for gen_paralog_node in
                                      node_cladogram.children if
                                      gen_paralog_node.is_leaf() == False):
                        if test_node.state == "y":
                            test_node_cptr += 1
                    if test_node_cptr == len(node_cladogram.children):
                        node_cladogram.state = "y"

                elif node_cladogram.clade_type == "ortholog":
                    orth_leaf_gen = (child_orth_leaf for child_orth_leaf in
                                     node_cladogram.children if
                                     child_orth_leaf.is_leaf() == True)
                    orth_leaf_cptr = len([orth_leaf for orth_leaf in
                                          orth_leaf_gen if
                                          orth_leaf.state == "y"])
                    orth_paralog_gen = (node_child_orth for node_child_orth in
                                        node_cladogram.get_descendants() if
                                    node_child_orth.clade_type == "paralog")
                    test_n = len([paralog_node for paralog_node in
                                  orth_paralog_gen if
                                  paralog_node.state == "n"])
                    if test_n == 0 and orth_leaf_cptr > 0:
                        node_cladogram.state = "y"

            for orth_node in (node for node in cladogram.traverse(
                    strategy="postorder") if node.clade_type =="ortholog"
                    and node.state == "y"):  # fill ortholog lineages with y
                for orth_node_child in orth_node.get_descendants():
                    orth_node_child.state = "y"

            # gather nodes in synapomorphies_set
            for leaf_syn_y in (leaf_syn for leaf_syn in cladogram.iter_leaves()
                               if leaf_syn.state == "y"):
                synapomorphies_set.add(find_synapomorphy(leaf_syn_y))

        # check nb of retained nodes: 1 (syn, accepted), 2 or more = homoplasy
        if len(synapomorphies_set) > 1:
            synapo_test = str("rejected")
        results_test_dict[char_names][char_states_names] = [synapo_test,
                                                            synapomorphies_set]

        return synapo_test, synapomorphies_set, results_test_dict, cladogram

    pdfs = []

    if rep_detector(cladogram_dict):
        print("ERROR: Repeated leaves have been detected in the " +
                         "cladogram.\nOperation aborted.")
        sys.exit(1)
    elif rep_detector(character_dict):
        print("ERROR: Repeated leaves have been detected in the " +
                         "character tree set.\nOperation aborted.")
        sys.exit(1)
    else:
        for t, i in character_dict.items():
            for leaf in t.get_leaf_names():
                if not leaf in list(cladogram_dict.keys())[0].get_leaf_names():
                    print("ERROR: Character tree " + str(i)
                                   + ": leaf set must be equal"
                                   + " or a subset of the cladogram leaf "
                                   + "set.\nOperation aborted.")
                    sys.exit(1)

    # activate stylenodes from ete3 if pdf option
    if pdf_files:
        try:
            from ete3 import NodeStyle, TreeStyle, faces, TextFace

            # node styles(standard, synapomorphy, homoplasy)
            nstyle_other = NodeStyle()  # neutral style
            nstyle_other["shape"] = "circle"
            nstyle_other["size"] = 0.3
            nstyle_other["fgcolor"] = "black"

            nstyle_accept = NodeStyle()  # synapomorphy style
            nstyle_accept["shape"] = "sphere"
            nstyle_accept["size"] = 3
            nstyle_accept["fgcolor"] = "limegreen"
            nstyle_accept.opacity = 0.3

            nstyle_reject = NodeStyle()  # homoplasy style
            nstyle_reject["shape"] = "sphere"
            nstyle_reject["size"] = 3
            nstyle_reject["fgcolor"] = "crimson"
            nstyle_reject.opacity = 0.3

            # leaf styles added dynamically (layout attached to the tree style)
            def lstyle(node_style):
                """
                Assigns layouts to nodes.

                """

                # assign leaf style
                if node_style.is_leaf():
                    if node_style.taxa_in_state == "in":
                        name_face = TextFace(node_style.name,
                                             fgcolor="steelblue",
                                             fsize=2,
                                             bold=True)
                    elif node_style.taxa_in_state == "out":
                        name_face = TextFace(node_style.name,
                                             fgcolor="dimgray",
                                             fsize=2)
                    else:
                        name_face = TextFace(node_style.name + " (?)",
                                             fgcolor="lightgray",
                                             fsize=2)
                    name_face.margin_left = 2
                    faces.add_face_to_node(name_face,
                                           node_style,
                                           column=0,
                                           position='branch-right')

        except ImportError:  # issue with ete3 imports
            print("A manual instal of PyQt5 is requested to use the --pdf "
                  + "functionality\nPlease install using 'pip install" +
                  " PyQt5' and rerun the command line")
            sys.exit(1)

    print("Initiating character state test procedure")

    # each node is numbered and defined as paralog/ortholog/apical
    for cladogram in cladogram_dict:
        cladogram, syn_dict = cladogram_label(cladogram,
                                              clade_number_option="yes",
                                              clade_type_option="yes",
                                              pdf_files=pdf_files)

    # list character, states, and taxa content + character cardinal
    character_cardinal_dict = {}
    character_component_dict = {}
    
    # for each character
    for character, values in character_dict.items():  
        character_states = {}
        character_states_count = 0
        value = str(values).split('.')[0] # remove artificial number from poly.
        
        # for each character state
        for node in character.traverse(strategy="preorder"):
            if node.is_leaf() == False  and node.is_root() == False:
                
                checkstate = True
                
                #find charstate names from hmatrix, new numbers if not present
                try: 
                    character_states_count = node.charstate_name  # hmatrix
                except AttributeError:
                    if type(character_states_count) == int:
                        character_states_count += 1
                        
                    # case of states without names in hmatrix
                    elif type(character_states_count) == str:
                        childlist = []
                        for children in node.children:
                            try:
                                childlist.append(children.charstate_name)
                            except:
                                pass
                        
                        character_states_count = "Unnamed?"
                        character_states_count += '+'.join(childlist)
                        
                        character_states["Character state #" + str(value) + "." 
                     + str(character_states_count)] = Tree.get_leaf_names(node)
                        
                        character_cardinal_dict["Character state #" 
                     + str(value) + "." + str(character_states_count)] = set(
                         Tree.get_leaf_names(character))
                        
                        checkstate = False
                        
                    else:
                        print("ERROR: incorrect label for character states")
                        sys.exit(1)
                
                if checkstate:
                
                    character_states["Character state #" + str(value) + "." 
                     + str(character_states_count)] = Tree.get_leaf_names(node)
                    
                    character_cardinal_dict["Character state #" + str(value) 
                                             + "." + str(
                                             character_states_count)] = set(
                                             Tree.get_leaf_names(character))
                
        character_component_dict["Character_"+str(value)  + "." + str(
                                 character_states_count)] = character_states

    # character state test and adding node style
    results_test_dict = {}
    for cladogram in cladogram_dict:
        with open(prefix+".chartest", "w") as results_file:
            results_file.write("#Character states tests")
            results_file.write(("\n#Legend: #character.state / state accepted"
                               " or rejected / node(s) characterised by the"
                               " state (if one: synapomorphy)"))
            results_file.write("\n" + cladogram.write(format=8) + "\n")

        # for each character
        for char_names, values in character_component_dict.items():
            results_test_dict[char_names] = {}

            #for each character state
            for char_states_names, taxa_in_state in values.items():

                # character state test
                if character_cardinal_dict[char_states_names] == set(
                        Tree.get_leaf_names(cladogram)):  # if no missing data
                    synapomorphy_test, synapomorphies_set, results_test_dict, \
                    cladogram = synapomorphy_test_placement(
                        cladogram, taxa_in_state, results_test_dict,
                        char_names, char_states_names)

                else:  # if there is missing data
                    pruned_cladogram = cladogram.copy(method='cpickle')
                    pruned_cladogram.prune(list(
                        character_cardinal_dict[char_states_names]))  #prune
                    pruned_cladogram, temp = cladogram_label(pruned_cladogram,
                             clade_number_option="no", clade_type_option="yes",
                             pdf_files=pdf_files)
                    synapomorphy_test, synapomorphies_set, results_test_dict, \
                        pruned_cladogram = synapomorphy_test_placement(
                            pruned_cladogram, taxa_in_state, results_test_dict,
                            char_names, char_states_names)

                    syn_list_extension = {}  # dict filled with synapomorphies
                    for synapomorphy in synapomorphies_set:
                        if synapomorphy.is_leaf():
                            syn_list_extension[synapomorphy] = [
                                synapomorphy.name]
                        elif synapomorphy.clade_type == "apical_node":
                            syn_list_extension[synapomorphy] = [leaf.name for
                            leaf in synapomorphy.iter_leaves() if leaf.name in
                            taxa_in_state]
                        else:
                            syn_list_extension[synapomorphy] = [leaf.name for
                            leaf in synapomorphy.iter_leaves()] #syn/leaf_list

                    synapomorphies_set_2 = set()
                    synapomorphies_set = set()

                    #dissociate internal nodes and leaves
                    for syn_list in syn_list_extension.values():
                        if len(syn_list) == 1:  # if leaf
                            synapomorphies_set_2.add(
                                cladogram.get_leaves_by_name(syn_list[0])[0])
                        else:
                            synapomorphies_set.add(
                                cladogram.get_common_ancestor(
                                    [cladogram.get_leaves_by_name(leaf_name)[0]
                                     for leaf_name in syn_list]))

                    # correction of y-leaves branched to apical nodes
                    for node in synapomorphies_set_2:
                        if node.get_ancestors()[0].clade_type == "apical_node":
                            synapomorphies_set.add(node.get_ancestors()[0])
                        elif (node.get_ancestors()[0].clade_type ==
                              "ortholog") and set([node_name.name for node_name
                              in node.get_ancestors()[0].get_leaves()]) & set(
                                  taxa_in_state) == {node.name}:
                            synapomorphies_set.add(node.get_ancestors()[0])
                        else:
                            synapomorphies_set.add(node)

                    # synapomorphies_set updated with the non pruned cladogram
                    results_test_dict[char_names][char_states_names] = [
                        synapomorphy_test, synapomorphies_set]

                # add leaf style for pdf files
                if pdf_files:

                    synapomorphy_style = TreeStyle()
                    synapomorphy_style.show_leaf_name = True
                    synapomorphy_style.margin_top = 10
                    synapomorphy_style.margin_right = 10
                    synapomorphy_style.margin_left = 10
                    synapomorphy_style.margin_bottom = 10
                    synapomorphy_style.show_scale = False
                    synapomorphy_style.show_leaf_name = False
                    synapomorphy_style.title.add_face(TextFace(
                        char_states_names, fsize=3), column=0) # pdf title

                    for node1 in cladogram.iter_leaves():
                        if node1.name in taxa_in_state:
                            node1.add_feature("taxa_in_state", "in")
                        elif node1.name in character_cardinal_dict[
                                char_states_names]:
                            node1.add_feature("taxa_in_state", "out")
                        else:
                            node1.add_feature("taxa_in_state", "?")

                    synapomorphy_style.layout_fn = lstyle # layout modif nodes

                    #ajout des styles de noeuds internes
                    for node_style in cladogram.traverse(strategy="postorder"):
                        if node_style in synapomorphies_set:
                            if synapomorphy_test == "accepted":
                                node_style.set_style(nstyle_accept)
                            else:
                                node_style.set_style(nstyle_reject)
                        else:
                            node_style.set_style(nstyle_other)

                    # add comments to the pdf
                    synapomorphy_style.legend.add_face(TextFace(
                        "Synapomorphy test: "+synapomorphy_test, fsize=2),
                        column=0)
                    synapomorphy_style.legend.add_face(TextFace(
                        "", fsize=2), column=0)
                    synapomorphy_style.legend.add_face(TextFace(
                        "Node(s): ", fsize=2), column=0)
                    node_line_pdf = str()
                    for node_text in sorted([syn_node.clade_label for syn_node
                                             in synapomorphies_set
                                             if not syn_node.is_leaf()]):
                        node_line_pdf += str(node_text)
                        node_line_pdf += " "
                    synapomorphy_style.legend.add_face(TextFace(
                        node_line_pdf, fsize=2), column=0)

                    node_line_pdf = str()
                    for node_text in sorted([syn_node.clade_label for
                                             syn_node in synapomorphies_set if
                                             syn_node.is_leaf()]):
                        node_line_pdf += str(node_text)
                        node_line_pdf += " "
                    synapomorphy_style.legend.add_face(TextFace(
                        node_line_pdf, fsize=2), column=0)

                    # save pdf
                    cladogram.render(char_states_names + ".pdf",
                                     tree_style=synapomorphy_style)

                    pdfs.append(char_states_names + ".pdf")

                # append text file
                with open(prefix+".chartest", "a") as results_file:
                    results_file.write("\n" + char_states_names + ": "
                                       + synapomorphy_test)
                    syn_not_leaf = sorted([syn_node.clade_label for syn_node in
                                           synapomorphies_set if not
                                           syn_node.is_leaf()])
                    syn_leaf = sorted([syn_node.clade_label for syn_node in
                                       synapomorphies_set if
                                       syn_node.is_leaf()])
                    if len(syn_not_leaf) == 1:
                        results_file.write(", clade ")
                    elif len(syn_not_leaf) > 1:
                        results_file.write(", clades ")
                    for node_syn in syn_not_leaf:
                        results_file.write(str(node_syn)+", ")
                    for node_syn in syn_leaf:
                        results_file.write(node_syn+", ")

    # write file with synapomorphies location
    for char_names, values1 in results_test_dict.items():
        for char_states_names, syn_test_set in values1.items():
            for node_set in syn_test_set[1]:
                if syn_test_set[0] == "accepted":
                    syn_dict[str(node_set.clade_label)]["accepted"].append(
                        [char_states_names.replace("Character state ","")])
                elif syn_test_set[0] == "rejected":
                    syn_dict[str(node_set.clade_label)]["rejected"].append([
                        char_states_names.replace("Character state ","")])

    with open(prefix+".chartest_node", "w") as results_file_tree:
        results_file_tree.write("Synapomorphies by node:")
        results_file_tree.write("\n" + cladogram.write(format=8) + "\n")
        for node_set in sorted((str(node_set) for node_set, values in
                                syn_dict.items())):
            if syn_dict[node_set]["accepted"] or (
                    syn_dict[node_set]["rejected"]):
                results_file_tree.write("\n# "+node_set)
                results_file_tree.write("\n     Synapomorphies : ")
                if syn_dict[node_set]["accepted"]:
                    for node_syn in syn_dict[node_set]["accepted"]:
                        results_file_tree.write(str(node_syn)+", ")
                results_file_tree.write("\n     Homoplasies : ")
                if syn_dict[node_set]["rejected"]:
                    for node_h in syn_dict[node_set]["rejected"]:
                        results_file_tree.write(str(node_h)+", ")
                results_file_tree.write("\n")

    if pdf_files:  # merge all pdfs

        merger = PdfMerger()

        for pdf in pdfs:
            merger.append(pdf)

        merger.write(prefix + ".pdf")
        merger.close()

        for pdf in pdfs:
            remove(pdf)

    print("Character state test procedure ended successfully")

    return results_test_dict, syn_dict


def describe_forest(character_dict, prefix, showtaxanames=False):
    """
    Write on a text file basic informations about rooted trees: number of
    nodes, terminals, internal nodes, symmetric nodes, apical nodes.
    Resolution of the tree, number of dichotomies, polytomies.
    Optionaly, the list of leaves can be writen on the file.

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
    prefix : str
        Prefix of the output file.
    showtaxanames : bool, optional
        If true, write the list of leaves of the tree in the output file.
        The default is False.

    Returns
    -------
    None.

    """

    def describe_tree(dtree, outfile, nb=1, showtaxanames=False):
        """
        Write basic informations about a single rooted tree.

        Parameters
        ----------
        dtree : ete3.Tree
            A single rooted tree.
        nb : integer
            number of the tree (default 1).

        Returns
        -------
        None.

        """

        terminals = sorted(dtree.get_leaf_names())
        cardinal = len(terminals)
        paralogn = 0
        orthologn = 0
        apicaln = 0
        totaln = 0
        dichotomies = 0
        polytomies = 0
        polytomiesd = dict()

        taxalist = taxa_extraction({dtree:0})
        repetitions = len(taxalist) - len(set(taxalist))

        for node in dtree.traverse("postorder"):
            if not node.is_leaf():  # if internal node

                totaln += 1
                i = len([child_node for child_node in node.get_children(
                    ) if not child_node.is_leaf()])
                j = len([child_node for child_node in node.get_children()])

                if i >= 2:  # paralog
                    paralogn += 1
                elif i == 1:  # ortholog
                    orthologn += 1
                elif i == 0:  # apical
                    apicaln += 1

                if j == 2:  # dichotomies vs polytomies
                    dichotomies += 1
                else:
                    polytomies += 1

                    if j in polytomiesd: # details about polytomies
                        polytomiesd[j] += 1
                    else:
                        polytomiesd[j] = 1

        with open(outfile, "a") as describefile:
            describefile.write("***************\n")
            describefile.write("Describe tree {}\n".format(str(nb)))
            describefile.write("Number of nodes: " + str(totaln + cardinal)
                               + "\n")
            describefile.write("    Terminals: " + str(cardinal)
                               + "\n")
            describefile.write("    Internal nodes: " + str(totaln)
                               + "\n")
            describefile.write("        Symmetric nodes: " + str(paralogn)
                               + "\n")
            describefile.write("        Asymmetric nodes: " + str(orthologn)
                               + "\n")
            describefile.write("        Apical nodes: " + str(apicaln)
                               + "\n")
            describefile.write("\n")

            if cardinal > 2:
                describefile.write("Resolution: " + str(
                    round(((totaln - 1) / (cardinal - 2)) * 100, 5)) + "%\n")
            else:
                describefile.write("Resolution: NA\n")

            describefile.write("Dichotomies: {}; polytomies: {}\n".format(
                str(dichotomies),str(polytomies)))
            describefile.write("Repetitions: {}\n".format(str(repetitions)))

            if polytomies > 0:
                describefile.write("Polytomies - details\n" +
                                   "(level: number of occurences)\n")
                for poly in sorted(polytomiesd.keys()):
                    describefile.write("{} : {}\n".format(str(poly),
                                                    str(polytomiesd[poly])))

            if showtaxanames:
                describefile.write("\nTaxa list:\n")
                for taxa in taxalist:
                    describefile.write(taxa+"\n")

            describefile.write("\n")

    i = 0
    for dtree in character_dict.keys():
        i += 1
        describe_tree(dtree, prefix+".dt", nb=i, showtaxanames=showtaxanames)

    print("Forest described")


def chartest(cladopath, charpaths, taxarep1=False, taxarep2=False, method="TMS", 
             prefix="AGATTA_chartest", pdf_files=False, verbose=False):
    """
    Run the character states testing procedure for hierarchical characters
    defined by Cao (2008) and corrected by Rineau (2017). This procedure tests
    at each symmetric node (paralog) of the cladogram the presence of
    instances of the derived state.

          Cao, N. (2008). Analyse à trois éléments et anatomie du bois des
          Fagales Engl (Doctoral dissertation, Paris, Muséum national
          d'histoire naturelle).

          Rineau, V. (2017). Un nouveau regard cladistique sur l'anatomie
          comparée, la phylogénie, la systématique et la paléoécologie des
          rudistes (Bivalvia, Hippuritida) (Doctoral dissertation,
          Université Pierre et Marie Curie).

    This function is made for giving trees directly to character_states_test
    from newick files

    Parameters
    ----------
    cladopath : str
        Path to a newick file containing the cladogram. The tree must be the
        optimal cladogram obtained from the cladistic analysis of character
        trees stored in charpaths.
    charpaths : list of str
        Path(s) to file(s) containing the initial character trees in newick or 
        hmatrix format.
    taxarep1 : str, optional
        DESCRIPTION. The default is False.
    taxarep2 : str, optional
        DESCRIPTION. The default is False.
    method : str, optional
        One of the two implemented algorithms of free-paralogy subtree
        analysis between "TMS" and "FPS" for removing polymorphism. 
        The default is "TMS".
    prefix : str
        Prefix of the files to be computed.
    pdf_files : str or bool, optional
        This argument can be the path stating the location where to save one
        pdf file for each character state for visualizing the results of the
        character state procedure.
        The default is False (no pdf files computed).
    verbose : bool, optional
        Verbose mode if True. The default is True.
        
    Returns
    -------
    None.

    """

    print("Loading cladogram")

    cladopath = path.expanduser(cladopath)
    charpaths = [path.expanduser(charpath) for charpath in charpaths]

    cladogram_dict = character_extraction(cladopath, taxarep1, verbose=False)

    print("Cladogram loaded")

    character_dict = wrapper_character_extraction(charpaths, 
                                                  taxarep2,
                                                  prefix,
                                                  verbose)
        
    # remove automatically repetitions if detected (user message printed)
    if rep_detector(character_dict):
        character_dict = del_replications_forest(character_dict,
                                                 method=method,
                                                 prefix=prefix,
                                                 verbose=verbose)

    # print(str(len(character_dict)) + " characters loaded")

    character_states_test(cladogram_dict, character_dict, prefix, pdf_files)
