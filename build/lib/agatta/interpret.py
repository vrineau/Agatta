# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:52 2020

@author: Valentin Rineau
"""

from fractions import Fraction
from itertools import combinations
from .ini import taxa_extraction
from .analysis import standard_tripdec
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
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

    constrict = treelist[0]  # base tree for the strict consensus computation
    del treelist[0]
    kill_node_list = []  # list of nodes to delete at the end

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
    Compute the reduced cladistic consensus (Wilkinson, 1994) from a list of
    rooted trees.

          Wilkinson, M. (1994). Common cladistic information and its consensus
          representation: reduced Adams and reduced cladistic consensus trees
          and profiles. Systematic Biology, 43(3), 343-368.

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
        raise ValueError("Trees must all have the same leaf set")

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

    #merging all nts into trees
    profile = []

    for cardinal, in_list in nts_groups.items():  #for each tree
        tree_profile = nts_to_trees(cardinal, in_list)
        if not len(tree_profile) == leaf_cptr:
            profile.append(tree_profile)

    profile.append(constrict(treelist, prefix=False, verbose=False))

    #sort trees by size
    def get_len(in_set):
        return len(in_set)

    profile = sorted(profile, key=get_len)

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


def RI(cladogram_dict, character_dict, weighting="FW", prefix=False):
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
        The tree must be the optimal cladogram obtained from the cladistic
        analysis of character trees stored in the character_dict argument.
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.
        The trees must be the initial characters.
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        . The default is "FW".
    prefix : str, optional
        Prefix of the saving file. The complete path can be
        used. The default is False (no file saved).

    Returns
    -------
    RI_char_dict : dict.
        Dictionary with character identifiers and global and per character RI.

    """

    print("Computing retention index")

    c_triplet_dict = standard_tripdec(cladogram_dict,
                                      weighting,
                                      prefix=False,
                                      verbose=False)
    RI_tot = [0, 0]
    RI_char_dict = {}
    RI_char_dict_num = {}
    RI_char_dict_denom = {}

    # computation of the global RI in the same time as RI per character
    for chartree, keys in character_dict.items():
        RI_char_dict[keys] = 0
        RI_char_dict_num[keys] = 0
        RI_char_dict_denom[keys] = 0
        triplet_dict = standard_tripdec({chartree: keys},
                                        weighting,
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

    for keys, values in RI_char_dict_num.items():
        if RI_char_dict_denom[keys] == 0:
            RI_char_dict[keys] = Fraction(0, 1)
        else:
            RI_char_dict[keys] = Fraction(RI_char_dict_num[keys],
                                          RI_char_dict_denom[keys])

    if RI_tot[1] == 0:
        RI_char_dict["Total retention index"] = Fraction(0, 1)  # add global RI
    else:
        RI_char_dict["Total retention index"] = Fraction(RI_tot[0], RI_tot[1])

    print("Retention index computed")

    # Write file
    if prefix:
        with open(prefix+".ri", "w") as RI_file:
            RI_file.write("#Retention index\n")

            for keys, values in RI_char_dict.items():  # IR list
                RI_string = str(keys)+": "+str(round((float(values)*100), 2))
                print(RI_string)
                RI_file.write("\n"+RI_string)

    return RI_char_dict


def ITRI(true_tree, reconstructed_tree, prefix, weighting="FW"):
    """
    Compute the inter-tree retention index (ITRI; Grand et al. 2014,
    which is an asymmetric distance between two trees, one
    reference tree and one tree to be compared to the reference tree.

    Grand, A., Zaragüeta-Bagils, R., Velez, L. M., & Ung, V. (2014).
    A cladistic re-analysis of the Gadiformes (Teleostei, Paracanthopterygii)
    using three-item analysis. Zootaxa, 3889(4), 525-552.

    Parameters
    ----------
    true_tree : dict
        Dictionary containing one newick tree (ete3 Tree object) as key.
        The tree is the reference tree.
    reconstructed_tree : dict
        Dictionary containing one newick tree (ete3 Tree object) as key to
        be compared to the reference tree.
    prefix : str
        Prefix of the text file to be saved containing the resulting ITRI.
        The complete path can be used. The default is False (no file saved).
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        The default is "FW".
    Returns
    -------
    float
        Resolving power in percentage.. Amount of information present in the
        compared tree present in the reference tree divided by the total
        information of the reference tree (Grand et al. 2014).
    float
        Artefactual resolution in percentage.. Amount of information present in
        the compared tree not present in the reference tree divided by the
        total information of the tree to be compared (Grand et al. 2014).
    float
        Efficiency. Resolving power - artefactual resolution in percentage.

    """

    tt_tripdic = standard_tripdec({true_tree: 0},
                                      weighting,
                                      prefix=False,
                                      verbose=False)

    rt_tripdic = standard_tripdec({reconstructed_tree: 0},
                                      weighting,
                                      prefix=False,
                                      verbose=False)


    RI_reconstructed_tree = Fraction(0, 1) #total reconstruted tree
    RI_true_tree = Fraction(0, 1) #total true tree
    RI_intersect_reconstructed_tree = Fraction(0, 1)
    RI_intersect_true_tree = Fraction(0, 1)

    for trip, FW in rt_tripdic.items():
        if trip in tt_tripdic:
            RI_intersect_reconstructed_tree += FW
        RI_reconstructed_tree += FW

    for trip, FW in tt_tripdic.items():
        if trip in rt_tripdic:
            RI_intersect_true_tree += FW
        RI_true_tree += FW

    ITRI_power = Fraction(RI_intersect_true_tree, RI_true_tree)

    try:
        ITRI_arte = 1 - Fraction(RI_intersect_reconstructed_tree,
                                 RI_reconstructed_tree)
    except ZeroDivisionError:
        ITRI_arte  = 0

    ITRI_efficiency = ((float(ITRI_power) * 100) - (float(ITRI_arte) * 100))

    if prefix:
        with open(prefix+".txt", "w") as itrifile:
            itrifile.write("power : " + str(round(float(ITRI_power), 3))
                           + "\n")
            itrifile.write("artefact : " + str(round(float(ITRI_arte), 3))
                           + "\n")
            itrifile.write("efficiency : "
                           + str(round(float(ITRI_efficiency), 3)) + "\n")

    print("power : " + str(round(float(ITRI_power), 3)) + "\n")
    print("artefact : " + str(round(float(ITRI_arte), 3)) + "\n")
    print("efficiency : " + str(round(float(ITRI_efficiency), 3)) + "\n")

    return float(ITRI_power)*100, float(ITRI_arte)*100, ITRI_efficiency


def triplet_distance(t1, t2, prefix,
                     method="itrisym_sum", weighting="FW"):
    """
    Compute the triplet distance between two trees. The triplet distance is the
    number of triplet differing between two trees (weights can be used).
    The order between the two trees t1 and t2 is not important.

    Parameters
    ----------
    t1 : dict
        Dictionary containing one newick tree (ete3 Tree object) as key.
        First tree.
    t2 : dict
        Dictionary containing one newick tree (ete3 Tree object) as key.
        Second tree.
    prefix : str
        Prefix of the text file to be saved containing the resulting triplet
        distance. The complete path can be used. The default is False
        (no file saved).
    method : TYPE, optional
        Two methods are available to compute the triplet distance:

            * 'itrisym_sum': (ITRI t1->t2 + ITRI t2->t1) / 2
            * 'itrisym_product': ITRI t1->t2 * ITRI t2->t1 (Grand et al. 2014)

        Grand, A., Zaragüeta-Bagils, R., Velez, L. M., & Ung, V. (2014).
        A cladistic re-analysis of the Gadiformes (Teleostei,
        Paracanthopterygii) using three-item analysis.
        Zootaxa, 3889(4), 525-552.

        The default is 'itrisym_sum'.
    weighting : str, optional
        Weighting scheme to use between:
                                 * FW (Fractional weighting from
                                       Rineau et al. 2021)
                                 * FWNL (Fractional weighting from
                                       Nelson and Ladiges),
                                 * MW (Minimal Weighting),
                                 * AW (Additive Weighting),
                                 * NW (No Weighting).
        The default is "FW".
    Returns
    -------
    ITRIsym : float
        Distance between the two trees.

    """

    t1_tripdic = standard_tripdec({t1: 0},
                                      weighting,
                                      prefix=False,
                                      verbose=False)
    t2_tripdic = standard_tripdec({t2: 0},
                                      weighting,
                                      prefix=False,
                                      verbose=False)

    RI_t2 = Fraction(0, 1)
    RI_t1 = Fraction(0, 1)
    RI_intersect_t2 = Fraction(0, 1)
    RI_intersect_t1 = Fraction(0, 1)

    for trip, FW in t2_tripdic.items():
        if trip in t1_tripdic:
            RI_intersect_t2 += FW
        RI_t2 += FW

    for trip, FW in t1_tripdic.items():
        if trip in t2_tripdic:
            RI_intersect_t1 += FW
        RI_t1 += FW

    ITRIp_12 = Fraction(RI_intersect_t1, RI_t1)
    ITRIp_21 = Fraction(RI_t1, RI_intersect_t1)

    try:
        ITRIa_12 = 1 - Fraction(RI_intersect_t2, RI_t2)
    except ZeroDivisionError:
        ITRIa_12  = 0

    try:
        ITRIa_21 = 1 - Fraction(RI_t2, RI_intersect_t2)
    except ZeroDivisionError:
        ITRIa_21  = 0

    ITRIe_12 = ((float(ITRIp_12)*100) - (float(ITRIa_12)*100))
    ITRIe_21 = ((float(ITRIp_21)*100) - (float(ITRIa_21)*100))

    if method == "itrisym_sum":  # classic mean
        ITRIsym = (ITRIe_12+ITRIe_21)/2

    elif method == "itrisym_product":  # Grand et al. 2014 proposal
        ITRIsym = ITRIe_12*ITRIe_21

    if prefix:
        with open(prefix+".txt", "w") as itrifile:
            itrifile.write(str(round(ITRIsym,3)))

    print("Symmetrical ITRI value :"+str(round(ITRIsym,3)))

    return ITRIsym


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
    pdf_files : str or bool, optional
        This argument can be the path stating the location where to save one
        pdf file for each character state for visualizing the results of the
        character state procedure.
        The default is False (no pdf files computed).

    Returns
    -------
    dict
        Two dictionaries with the results of the procedure sorted in two ways.

    """

    # label cladogram nodes (types and id)
    def cladogram_label(cladogram, clade_number_option, clade_type_option):
        """
        Label tree nodes.
        """

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

        except ImportError as error:  # issue with ete3 imports
            print("""Installing PyQt5 is requested to use this functionality\n
                  To fix the issue, please read :
                  https://github.com/etetoolkit/ete/issues/354""" + str(error))

    print("Initiating character state test procedure")

    # each node is numbered and defined as paralog/ortholog/apical
    for cladogram in cladogram_dict:
        cladogram, syn_dict = cladogram_label(cladogram,
                                              clade_number_option="yes",
                                              clade_type_option="yes")

    # list character, states, and taxa content
    character_component_dict = {}
    for character, values in character_dict.items():
        character_states = {}
        character_states_count = 1
        for node in character.traverse(strategy="preorder"):
            if node.is_leaf() == False  and node.is_root() == False:
                character_states["Character state #" + str(values) + "." + str(
                           character_states_count)] = Tree.get_leaf_names(node)
                character_states_count += 1
        character_component_dict["Character_"+str(values)] = character_states
    # character cardinal
    character_cardinal_dict = {}
    for character, values in character_dict.items():
        character_states_count = 1
        for node in character.traverse(strategy="preorder"):
            if node.is_leaf() == False  and node.is_root() == False:
                character_cardinal_dict["Character state #" + str(values) + "."
                                        + str(character_states_count)] = set(
                                            Tree.get_leaf_names(character))
                character_states_count += 1

    # character state test and adding node style
    results_test_dict = {}
    for cladogram in cladogram_dict:
        with open(prefix+".chartest", "a") as results_file:
            results_file.write("#Character states tests")
            results_file.write("""\n#Legend: #character.state / state accepted
                               or rejected / node(s) characterized by the state
                               (if one: synapomorphy)""")
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
                             clade_number_option="no", clade_type_option="yes")
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
                    cladogram.render(pdf_files+char_states_names + ".pdf",
                                     tree_style=synapomorphy_style)


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
            describefile.write("Repetitions:{}\n".format(str(repetitions)))

            if polytomies > 0:
                describefile.write("""Polytomies - details
                                   (level: number of occurences)\n""")
                for poly in sorted(polytomiesd.keys()):
                    describefile.write("{} : {}\n".format(str(poly),
                                                    str(polytomiesd[poly])))

            if showtaxanames:
                describefile.write("\nTaxa list:\n")
                for taxa in taxalist:
                    describefile.write(taxa+"\n")


    i = 0
    for dtree in character_dict.keys():
        i += 1
        describe_tree(dtree, prefix+".dt", nb=i, showtaxanames=showtaxanames)
