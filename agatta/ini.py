# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:46:58 2020

@author: Valentin Rineau
"""

from os import path
from re import search
from re import findall
from re import compile
from tkinter import Tk
from tkinter import filedialog
from collections import defaultdict
import warnings
import csv

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree


def character_extraction(infile=False, taxa_replacement=False):
    """
    Extract newick rooted trees from a text file.
    The function extract lines begining with an open parenthesis and
    closing with a semicolon and interpret them as newick strings.
    The first argument is the path of the file, but the function can be used
    without any argument. In this case, a file selection window appears.
    A taxa bloc can also be present to replace taxa names.

    Parameters
    ----------
    infile : str, optional
        Path of the file containing a single newick tree on each line.
        Example:

            (a,b,(c,d,(e,f)));
            ((c,e),(a,(b,(d,f))));

        The default is False (a selection window appears in this case).

    taxa_replacement : bool, optional
        If set to True, will replace taxa according to a taxa bloc, i.e.,
        set of lines in the form new_name = old_name. The separator
        between names must be ' = '. Example:
             AA = Diceras
             AB = Valletia
             AC = Monopleura
        The default is False.

    Returns
    -------
    character_dict : dict
        Dictionnary containing trees in keys with indices in values.

    """

    print("Loading character trees")

    character_dict = {}  # return dictionary containing trees
    taxa_dict = {}  # taxa / symbol dictionary used if taxa_replacement == True

    # selection file window
    if infile is None:
        root = Tk()
        infile = filedialog.askopenfilename(
            title="Select file containing newick character trees",
            filetypes=(
                ("all files", "*.*"),
                ("tree files", "*.tre"),
                ("3ia Lisbeth 1.0 files", "*.3ia")))
        root.withdraw()

    fileName, fileExtension = path.splitext(infile)

    # open input file and extract newick strings
    with open(infile, "r") as file_tree:
        a = 1
        line_nb = 0
        for line in file_tree:
            line_nb += 1
            for character_newick in findall(r"\([^ \t\n\r\f\v]+;", line):
                if character_newick:

                    try:
                        character_dict[Tree(character_newick)] = a
                        a += 1

                    except Exception as error:
                        raise UserWarning("Error in the file " + infile +
                                          "\nLine " + str(line_nb) + ": " +
                                          str(error))

            # searching taxa bloc for replacing taxa names
            if taxa_replacement and search(r"\s=\s", line):
                taxa_dict[line.split(" = ")[0].split()[0]] = line.split(
                                                            " = ")[1].strip()

    if taxa_dict:
        for cladogram in character_dict.keys():
            for leaf in cladogram.iter_leaves():
                if taxa_dict[leaf.name]:
                    leaf.name = taxa_dict[leaf.name]

    print("{} characters loaded".format(str(len(character_dict))))

    return character_dict

def taxa_extraction(character_dict):
    """
    Extracts all leaves from a dictionary containing newick trees in keys.

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.

    Returns
    -------
    taxa_set : set
        set of taxa names.

    """

    taxa_set = set()
    for character in character_dict.keys():
        taxa_set.update(set(Tree.get_leaf_names(character)))
    return taxa_set


def taxa_triplet_extraction(triplet_dict):
    """
    Extracts all leaves from a dictionary containing triplets in keys.

    Parameters
    ----------
    triplet_dict : dict
        Dictionary containing triplets (agatta triplet objects) in keys.

    Returns
    -------
    taxa_set : set
        set of taxa names.

    """

    taxa_set = set()

    for triplet in triplet_dict.keys():
        taxa_set.update(triplet.out_taxa)
        taxa_set.update(triplet.in_taxa)

    return taxa_set


def taxa_to_numbers(character_dict):
    """
    Extracts all leaves from a dictionary containing newick trees in keys,
    and assign to each taxa name an identifier.

    Parameters
    ----------
    character_dict : dict
        Dictionary containing newick trees (ete3 Tree objects) in keys.

    Returns
    -------
    taxa_dict : dict
        Dictionary with taxa names as keys and identifiers as values.
    code_dict : dict
        Dictionary inverse of taxa_dict used for translation of identifiers.
    taxa_convers : list
        list of strings with taxa names and corresponding id.

    """
    taxa_set = taxa_extraction(character_dict)
    taxa_dict = {}
    code_dict = {}
    taxa_convers = []
    taxa_cptr = 1

    for taxa1 in taxa_set:
        taxa_dict[taxa1] = taxa_cptr
        code_dict[taxa_cptr] = taxa1
        taxa_convers.append(str(taxa1)+" "+str(taxa_cptr))
        taxa_cptr += 1

    return taxa_dict, code_dict, taxa_convers


def standardisation(tree_file, biogeo_tab, prefix, verbose=False):
    """
    This function is used for standardisation of characters in the context of
    cladistic biogeography, i.e., replacement of leaves using a correspondence
    table. In cladistic theory, the standardisation is the construction of
    character trees (hypotheses of kinship relationship between bearers) based
    on homologies (hypotheses of kinship relationships between parts). The
    standardisation is used to convert phylogenies in areagrams in cladistic
    biogeography.
    MAST are automatically removed.
    A file named prefix.std is generated at the end.

    Parameters
    ----------
    tree_file : str, optional
        Path of the file containing a single newick tree on each line.
        Example:

            (a,b,(c,d,(e,f)));
            ((c,e),(a,(b,(d,f))));

        The default is False (a selection window appears in this case).
    biogeo_tab : str
        path of table with semicolon separators and two columns, taxa in
        the left column and corresponding areas in the right column.
    prefix : str, optional
        prefix of the .std file to save.
    verbose : bool, optional
        If true, print in the terminal the resulting areagrams and informations
        on taxa, areas, MASTs (Multiple Area Single Taxa) and repetitions.
        The default is False.

    Returns
    -------
    areagram_dict : dict
        Dictionary with areagram trees in keys (ete Tree objects) and same
        index as character_dict to track the original trees.

    """

    def biogeo_table(biogeo_tab, verbose=False):
        """
        This function requires in input a table (semicolon separators) with
        the correspondence between parts in the left column (e.g., taxa in
        cladistic biogeography) and bearers in the right column (e.g., areas)

        Parameters
        ----------
        biogeo_tab : str
            path of table with semicolon separators and two columns, taxa in
            the left column and corresponding areas in the right column.
        verbose : bool, optional
            If true, print in the terminal informations on taxa, areas, MASTs
            (Multiple Area Single Taxa) and repetitions. The default is False.

        Returns
        -------
        biogeo_dict : dict
            dictionary with taxa in keys and areas in values.
        MAST : set
            set of MAST areas.
        taxrep : set
            set of repeated areas.

        """

        biogeo_dict = {}
        taxa = set()    # taxa set
        areas = set()   # areas set
        MAST = set()    # set of detected MAST
        taxrep = set()  # set repeated leaves

        with open(biogeo_tab, "r") as bt_file:
            for line in bt_file:
                regex = compile(r'[\n\r\t]')
                line = regex.sub("", line)
                linesplt = line.split(";")

                # detection of MASTs (if taxon already in list)
                if linesplt[0] in taxa:
                    biogeo_dict[linesplt[0]].append(linesplt[1])
                    MAST.add(linesplt[0])

                else:
                    biogeo_dict[linesplt[0]] = [linesplt[1]]

                # detection of repeated areas
                if linesplt[1] in areas:
                    taxrep.add(linesplt[1])

                taxa.add(linesplt[0])
                areas.add(linesplt[1])

        if verbose:
            print(str(len(taxa)) + " terminal taxa, " + str(len(areas)) +
                  " terminal areas")
            print(str(len(taxrep)) + " repeated taxa, " + str(len(MAST)) +
                  " MAST detected")
            print("taxa: {}".format(taxa))
            print("areas: {}".format(areas))
            if len(taxrep) > 0:
                print("repeated taxa: {}".format(taxrep))
            if len(MAST) > 0:
                print("MAST: {}".format(MAST))

        return biogeo_dict, MAST, taxrep

    biogeo_dict, MAST, taxrep = biogeo_table(biogeo_tab, verbose=False)
    if type(tree_file) == str:
        character_dict = character_extraction(tree_file)
    else:
        character_dict = tree_file
    areagram_dict = dict()
    nb_char = len(character_dict)

    if nb_char > 1:
        print("Standardising {} characters".format(str(nb_char)))
    else:
        print("Standardising 1 character")

    # for each tree
    with open(prefix + ".std", "w") as area_file:
        for phylogeny, index in character_dict.items():

            phylogeny2 = phylogeny.copy()

            # MAST deletion
            for leaf in phylogeny.get_leaves():
                if leaf.name in MAST:
                    leaveswmast = list(set(phylogeny2.get_leaf_names()) - MAST)
                    phylogeny2.prune(leaveswmast)

            # for each terminal taxa, replacement by its corresponding area
            # MAST excepted
            for leaf in phylogeny.get_leaves():
                if leaf.name not in MAST:
                    for leafnode in phylogeny2.get_leaves_by_name(leaf.name):
                        leafnode.name = biogeo_dict[leaf.name][0]

            # add new areagram
            areagram_dict[phylogeny2] = index
            new_tree_line = str(index)
            new_tree_line += "    "
            new_tree_line += phylogeny2.write(format=9)
            new_tree_line += "    "
            new_tree_line += phylogeny.write(format=9)

            area_file.write(new_tree_line)

            if verbose:
                print(new_tree_line)

    print("Characters standardised")

    return areagram_dict


def hmatrix(path, prefix=False, chardec=False, verbose=False):
    """
    Function that build ete3 trees from a file containing a hierarchical
    matrix. The matrix must be a table text (e.g., a csv file) with semicolons
    as separators.

    Format of the matrix:

      * The first line corresponds to the hierarchical structure of the
        character tree with labels.
        It is a newick tree without ending semicolon and without leaves.
        It is to the user to choose which nodes are labelled. Labels must be
        integers.

        Examples: (0(1)), (0(1)(2(3(4)))), ((0)(1)).

      * The first column corresponds to taxa names. Semicolons cannot be used
        in a name. The first cell is empty as it corresponds to the line
        of hierarchical characters.

      * The following lines are identical to those of a standard
        taxon/character matrix. Each state corresponds to the branching
        position of a specific taxon on the character.
        Each cell can be filled by:

            - A question mark "?" if the state for a taxon is unknown.

            - An integer which should be mandatory present in the first line
              of the matrix. For example, if the character is coded as (0,(1))
              by the user, only 0 and 1 can be used.

            - In case of polymorphism (if the user want to branch several
              instances of the same taxon in several nodes), the several
              integers must be separated by a comma ",".

        Example:

            ;(0,(1));(0,(1));(0,(1,(2)));(0,(1,(2),(3)))
            Diceras;0;0;0;0
            Hippurites;0;0;0;1
            Radiolites;1;0;2;2
            Titanosarcolites;1;1;2;2,3
            Clinocaprina;1;1;1;3
            Vaccinites;0,1;1;2;3


    Parameters
    ----------
    path : str
        path of the file containing the hierarchical matrix.
    prefix : str, optional
        Prefix of the file to save. If prefix is set to false, no file is
        saved. The default is False.
    chardec : bool, optional
        If true, decompose each tree into components (one tree for each
        informative node). The default is False.
    verbose : bool, optional
        If true, print all generated trees. The default is False.

    Returns
    -------
    character_dict : dict
        Dictionnary containing trees in keys with indices in values.

    """

    print("Loading hierarchical matrix")

    character_dict = dict()  # trees without polytomies
    temp_character_dict = dict()  # trees with polytomies (raw data)
    error_dict = defaultdict(list)

    # read matrix
    with open(path, 'rt') as f:
        data = csv.reader(f, delimiter=';')
        hmatrix = list(data)

    print("Hierarchical matrix loaded")
    print("Treefication of the hierarchical matrix.")

    # construction of character trees backbone
    i = 1
    treeliststr = hmatrix[0]
    del treeliststr[0]

    for char in treeliststr:
        temp_character_dict[Tree(char+";")] = str(i)
        i += 1

    taxalist = [hmatrix[ncol][0] for ncol in range(1, len(hmatrix))]

    # option for character decomposition into components
    if chardec:
        temp_character_dict_binary = dict()

        for char, ind in temp_character_dict.items():  # for each character

            # for each unrooted character state
            for state in char.traverse(strategy="levelorder"):
                if state.is_leaf() and not state.get_ancestors()[0].is_root():

                    statetree = char.copy()

                    for isolstate in statetree.traverse(strategy="levelorder"):
                        if not isolstate.is_leaf() and not isolstate.is_root():
                            if not isolstate.get_children()[0].name == \
                                                                    state.name:
                                isolstate.delete()

                    stateline = str(ind) + "." + str(state)
                    temp_character_dict_binary[statetree] = stateline

        temp_character_dict = temp_character_dict_binary

    # fill characters
    for char, ind in temp_character_dict.items():  # for each character

        # mark branches to delete at the end
        delnodes = []  # list of leaves to delete
        for leaf in char.iter_leaves():
            delnodes.append(leaf)
            leaf.name = "_agatta_charstate_"+leaf.name  # for avoiding the
            # specific case where character state and taxa have the same name

        # for each taxon
        for taxa in taxalist:
            if chardec:
                ind2 = ind.split(".")[0]
            else:
                ind2 = ind
            charstates = hmatrix[taxalist.index(taxa)+1][int(ind2)].split(",")

            if charstates[0] != "?":  # if not missing data

                # branching
                for charstate in charstates:

                    # if incorrect matrix, raise error, otherwise, connect taxa
                    try:
                        branchnode = char.get_leaves_by_name(
                            "_agatta_charstate_" +
                            charstate)[0].get_ancestors()[0]
                        branchnode.add_child(name=taxa)
                    except:
                        error_dict[ind2].append(taxa)  # id character / taxon

        # deletion of state branches
        for delnode in delnodes:
            delnode.delete()

    for char, value in temp_character_dict.items():
        char.ladderize()
        character_dict[char] = str(value)

    # mode verbose a transformer en raise error
    if verbose:
        if error_dict:

            print("Error in the input matrix\n")

            for charnum, taxalist in error_dict.items():
                print("Character number {}:\n".format(str(charnum)))
                for taxa in taxalist:
                    print(" - "+taxa)

    # save resulting tree file
    if prefix:
        with open(prefix+".hmatrix", "w") as treefile:
            for char, charnum in character_dict.items():
                treefile.write(str(charnum)+" : "+char.write(format=9)+"\n")

    print("{} characters computed from the matrix".format(
                                                    str(len(character_dict))))

    return character_dict


def helper(command):
    """
    Helper function that print the help block corresponding to a command.

    Parameters
    ----------
    option : str
        Any command that can be used with agatta. Agatta commands comprises
        analysis, tripdec, ri, chartest, convert, fp, consensus, describetree,
        standardisation, hmatrix

    Returns
    -------
    None.

    """
    if command == "analysis":
        print("Help block line for analysis command")

    elif command == "tripdec":
        print("Help block line for tripdec command")

    elif command == "ri":
        print("Help block line for ri command")

    elif command == "chartest":
        print("Help block line for chartest command")

    elif command == "convert":
        print("Help block line for convert command")

    elif command == "fp":
        print("Help block line for fp command")

    elif command == "consensus":
        print("Help block line for consensus command")

    elif command == "describetree":
        print("Help block line for describetree command")

    elif command == "standardisation":
        print("Help block line for standardisation command")

    elif command == "hmatrix":
        print("Help block line for hmatrix command")
