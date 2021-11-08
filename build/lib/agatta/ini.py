# -*- coding: utf-8 -*-
"""

    AGATTA: Three-item analysis Python package
    By Valentin Rineau and Paul Zaharias

    AGATTA is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from os import path
from re import search
from re import findall
from tkinter import Tk
from tkinter import filedialog
import csv
import sys
import warnings
import treeswift

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    from ete3 import Tree
    from ete3.parser.newick import NewickError


def character_extraction(infile=False, taxa_replacement=False, verbose=True):
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
        Path of the file containing rooted trees in newick format
        (https://en.wikipedia.org/wiki/Newick_format). The input file must
        contain a single newick string on each line. Whitespaces are not
        allowed for leaf names.
        Example:

            (a,b,(c,d,(e,f)));
            ((c,e),(a,(b,(d,f))));

        The default is False (a selection window appears in this case).

    taxa_replacement : str, optional
        Path of a table file containing two columns. The first column
        corresponds to the names of the terminals of the newick stored in
        infile, and the second column corresponds to their names the user wants
        to obtain at the end. All separators accepted. Example:
             AA Diceras
             AB Valletia
             AC Monopleura
        The default is False (no replacement).

    Returns
    -------
    character_dict : dict
        Dictionnary containing trees in keys with indices in values.

    """

    if verbose:
        print("Loading character trees")

    character_dict = {}  # return dictionary containing trees
    taxa_dict = {}  # taxa / symbol dictionary used if taxa_replacement == True

    # selection file window
    if infile is None:
        root = Tk()
        infile = filedialog.askopenfilename(
            title="Select file containing newick character trees",
            filetypes=(
                ("Newick files", "*.*"),
                ("Lisbeth files", "*.3ia"),
                ("Nexus files", "*.nex"),
                ("Nexml files", "*.nexml")))
        root.withdraw()

    if not path.isfile(infile):
        sys.exit(print("ERROR: The file '" + infile + "' does not exist." +
                       "\nOperation aborted."))

    fileName, fileExtension = path.splitext(infile)
    a = 1

    if fileExtension == ".3ia":
        with open(infile, "r") as file_tree :
            line_nb = 0
            for line in file_tree:
                line_nb += 1
                for character_newick in findall(r"\(\S+;", line):
                    if character_newick:
                        try:
                            character_dict[Tree(character_newick)] = a
                            a += 1
                        except NewickError:
                            sys.exit(print("ERROR: " +
                                "Line {}: Broken newick structure.".format(
                                str(line_nb))))

                if search(r"\s=\s", line):
                    taxa_dict[line.split(" = ")[0].split()[0]] = line.split(
                        " = ")[1].strip()

    else:
        try:
            if fileExtension == ".nex":
                tstreelist = treeswift.read_tree_nexus(infile)

            elif fileExtension == ".nexml":
                tstreelist = treeswift.read_tree_nexml(infile)

            else:  # newick files
                tstreelist = treeswift.read_tree_newick(infile)
        except:
            sys.exit(print("ERROR: The file is broken. Please check the" +
                           "format. \nNexus files must have the .nex " +
                           "extension.\nNexml files must have the .nexml " +
                           "extension.\nThe Lisbeth input files must have " +
                           "the .3ia extension.\nAll other extension are " +
                           "considered as newick files containing only " +
                           "newick strings (one on each line)."))

        if not isinstance(tstreelist, list):
            tstreelist = [tstreelist]

        for tstree in tstreelist:
            if isinstance(tstree, dict):
                for idtree, tst in tstree.items():
                    try:
                        if (tst.newick().startswith("[&R] ") or
                            tst.newick().startswith("[&U] ")):
                            character_dict[Tree(
                                tst.newick().split(" ")[1])] = idtree
                        else:
                            character_dict[Tree(tst.newick())] = idtree
                    except NewickError:
                        sys.exit(print(
                            "Tree {}: Broken newick structure.".format(
                                str(idtree))))

            else:
                try:
                    if (tstree.newick().startswith("[&R] ") or
                        tstree.newick().startswith("[&U] ")):
                        character_dict[Tree(tstree.newick().split(" ")[1])] = a
                    else:
                        character_dict[Tree(tstree.newick())] = a
                    a += 1
                except NewickError:
                    sys.exit(print(
                        "Tree {}: Broken newick structure.".format(str(a))))

    if taxa_replacement:
        if not path.isfile(taxa_replacement):
            sys.exit(print("ERROR: The input file '" + taxa_replacement
                                    + "' does not exist."
                                    + "\nOperation aborted."))

        with open(taxa_replacement, "r") as taxa_table:
            try:
                dialect = csv.Sniffer().sniff(taxa_table.read())
            except:
                sys.exit(print("ERROR: Could not determine separator in the"
                               + " table file '" + taxa_replacement
                               + "'.\nThe table is probably broken."
                               + "\nOperation aborted."))

        with open(taxa_replacement, "r") as taxa_table:
            tab_test = list(csv.reader(taxa_table,
                                       delimiter=dialect.delimiter))
            if not len({len(l) for l in tab_test}) == 1:
                 sys.exit(print("ERROR: The table file '"
                                   + taxa_replacement + "' is broken."
                                   + "\nOperation aborted."))

            for idtax, nametax in tab_test:
                taxa_dict[idtax] = nametax

        for cladogram in character_dict.keys():
            for leaf in cladogram.iter_leaves():
                try:
                    leaf.name = taxa_dict[leaf.name]
                except KeyError:
                    sys.exit(print("ERROR: The name '" + str(leaf.name)
                                   + "' does not exists in the table file '"
                                   + taxa_replacement
                                   + "'.\nOperation aborted."))

    if verbose:
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
    tree_file : str
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
            path of table file with two columns, taxa in
            the left column and corresponding areas in the right column.
            Example:

                Clinocaprina USA
                Titanosarcolites Jamaica
                Diceras France

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
            try:
                dialect = csv.Sniffer().sniff(bt_file.read())
            except:
                sys.exit(print("ERROR: The table file '" + biogeo_tab +
                                  "' is probably broken. Could not determine" +
                                  "separator.\nOperation aborted."))

        with open(biogeo_tab, "r") as bt_file:

            tab_test = csv.reader(bt_file, delimiter=dialect.delimiter)
            table = list(tab_test)
            if not len({len(l) for l in table}) == 1:
                 sys.exit(print("ERROR: The table file '" + biogeo_tab +
                                  "' is broken.\nOperation aborted."))

        for line in table:
            # detection of MASTs (if taxon already in list)
            if line[0] in taxa:
                biogeo_dict[line[0]].append(line[1])
                MAST.add(line[0])

            else:
                biogeo_dict[line[0]] = [line[1]]

            # detection of repeated areas
            if line[1] in areas:
                taxrep.add(line[1])

            taxa.add(line[0])
            areas.add(line[1])

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
                        try:
                            leafnode.name = biogeo_dict[leaf.name][0]
                        except KeyError:
                            sys.exit(print("ERROR: The name '" +
                                             str(leaf.name) +
                                             "' does not exists in the table" +
                                             " file '" + biogeo_tab + "'.\n" +
                                             "\nOperation aborted."))

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


def hmatrix(infile, prefix=False, chardec=False, verbose=False):
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
    infile : str
        Path of the file containing the hierarchical matrix.
    prefix : str, optional
        Prefix of the file to save. If prefix is set to false, no file is
        saved. The default is False.
    chardec : bool, optional
        If true, decompose each tree into components (one tree for each
        informative node). The default is False.

    Returns
    -------
    character_dict : dict
        Dictionnary containing trees in keys with indices in values.

    """

    print("Loading hierarchical matrix")

    character_dict = dict()  # trees without polytomies
    temp_character_dict = dict()  # trees with polytomies (raw data)

    # read matrix
    if not path.isfile(infile):
        sys.exit(print("ERROR: The hierarchical matrix '" + infile
                                + "' does not exist.\nOperation aborted."))

    with open(infile, 'rt') as f:
        try:
            dialect = csv.Sniffer().sniff(f.read())
        except:
            sys.exit(print(sys.exit(print("ERROR: The table file '" + infile +
                                  "' is probably broken. Could not determine" +
                                  "separator.\nOperation aborted."))))

    with open(infile, 'rt') as f:
        if dialect.delimiter == ",":
            sys.exit(print("ERROR: ',' can't be a delimiter in the " +
                           "hierarchical matrix format(usage restricted for " +
                           "polymorphic instances).\nOperation aborted."""))

        if dialect.delimiter not in [",",";","\t"," ","|"]:
             sys.exit(print("""ERROR: Error in the hierarchical matrix format.
                               The separator must be one of these:

                                  - semicolon ';'
                                  - tabulation '   '
                                  - space ' '
                                  - pipe '|'

                                Operation aborted."""))

        data = csv.reader(f, delimiter=dialect.delimiter)
        hmatrix = list(data)
        if not len({len(l) for l in hmatrix}) == 1:
             sys.exit(print("ERROR: the hierarchical matrix '" +
                            infile + "' is broken.\nOperation aborted."))

    print("Hierarchical matrix loaded")
    print("Treefication of the hierarchical matrix.")

    # construction of character trees backbone
    i = 1
    treeliststr = hmatrix[0]
    del treeliststr[0]

    for char in treeliststr:
        try:
            temp_character_dict[Tree(char+";")] = str(i)
        except NewickError:
            sys.exit(print("ERROR: The newick tree " + char + " in column " +
                               str(i+1) + "is broken\nOperation aborted."))
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
                        sys.exit(print("ERROR: The instance '" + str(ind2)
                                         + "' of leaf " + taxa + " in column "
                                         + str(ind) + " does not match.\n"
                                         +"Operation aborted."))

        # deletion of state branches
        for delnode in delnodes:
            delnode.delete()

    for char, value in temp_character_dict.items():
        char.ladderize()
        character_dict[char] = str(value)

    # save resulting tree file
    if prefix:
        with open(prefix+".hmatrix", "w") as treefile:
            for char, charnum in character_dict.items():
                treefile.write(str(charnum)+" : "+char.write(format=9)+"\n")

    print("{} characters computed from the matrix".format(
                                                    str(len(character_dict))))

    return character_dict


def checkargs(arguments):
    """
    Check parsing and exit if error in flags

    Parameters
    ----------
    arguments : dict
        Docopt dict.

    Returns
    -------
    None.

    """
    if arguments["analysis"]:

        if arguments.get("--pdf", False) and not arguments["--chartest"]:
            sys.exit(print("ERROR: --pdf flag cannot be used without " +
                           "--chartest"))

        if (arguments["--software"] == "paup" or
            arguments["--software"] == "tnt"):
            if not arguments.get("--softpath", False):
                sys.exit(print("ERROR: For PAUP and TNT, the path of the " +
                      "software must be completed in --softpath"))

        elif arguments["--software"] != "agatta":
            sys.exit(print("ERROR: --software flag must be one of: 'paup', " +
                           "'tnt', 'agatta'"))

    if arguments["convert"]:

        if not arguments["--filetype"] in ["trees","triplets"]:
            sys.exit(print("ERROR: --filetype flag must be one of: 'trees', " +
                           "'triplets'"))

        try:
            int(arguments["--multiplier"])
        except ValueError:
            sys.exit(print("ERROR: --multiplier flag must be an integer"))

    if arguments["convert"] or arguments["analysis"]:

        if not arguments["--analysis"] in ["auto","heuristic", "bandb"]:
            sys.exit(print("ERROR: --analysis flag must be one of: 'auto', " +
                           "'heuristic', 'bandb"))

        try:
            int(arguments["--nrep"])
        except ValueError:
            sys.exit(print("ERROR: --nrep flag must be an integer (number of" +
                  " replicates in heuristic search"))

    if arguments["tripdec"] or arguments["convert"] or arguments["analysis"]:

        if arguments["--parallel"] != "auto":
            try:
                int(arguments["--parallel"])
            except ValueError:
                sys.exit(print("ERROR: --parallel flag must be 'auto' or an " +
                      "integer specifying the number of cpu"))

    if arguments["fp"] or arguments["analysis"]:
        if not arguments["--method"] in ["Rineau","Nelson"]:
            sys.exit(print("ERROR: --method flag must be one of: 'Rineau', " +
                           "'Nelson'"))

    if arguments["consensus"] or arguments["analysis"]:
        if arguments.get("--consensus", False):
            if not arguments["--consensus"] in ["strict","rcc"]:
                sys.exit(print("ERROR: --consensus flag must be one of: " +
                               "'strict', 'rcc'"))

    if arguments["support"]:
        if not arguments["--index"] in ["ri","itri","itrisym_sum",
                                        "itrisym_product"]:
            sys.exit(print("ERROR: --index flag must be one of: 'ri', " +
                               "'itri', 'itrisym_sum', 'itrisym_product"))


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

    else:
        print("""This command does not exist. The current commands are:
              analysis, tripdec, ri, chartest, convert, fp, consensus,
              describetree, standardisation, hmatrix""")