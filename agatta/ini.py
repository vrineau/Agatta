# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    By Valentin Rineau and Paul Zaharias

    Agatta is a set of tools in the cladistic framework to perform
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
        Dictionary containing triplets (Agatta triplet objects) in keys.

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
                treefile.write(str(charnum) + " : "
                               + char.write(format=9) + "\n")

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
            int(arguments["--replicates"])
        except ValueError:
            sys.exit(print("ERROR: --replicates flag must be an integer " +
                  "(number of replicates in heuristic search"))

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
        Any command that can be used with Agatta. Agatta commands comprises
        analysis, tripdec, ri, chartest, convert, fp, consensus, describetree,
        standardisation, hmatrix

    Returns
    -------
    None.

    """
    if command == "analysis":
        print("\033[1;34;48m    analysis")
        print("""
    Main function of the Agatta python package. Allow to perform a three-item
    analysis, e.g., in the context of systematics phylogenetics, using
    hierarchical characters, or in cladistic biogeography. The cladogram(s)
    obtained by congruence is the tree that maximises the amount of hypotheses
    of cladistic relationships (i.e., the three-item statements) deduced from
    the input trees (the characters).

    The analysis can be performed using a text file containing a hierarchical
    matrix (see section mandatory parameters below for informations about
    the format) or newick rooted trees encoded in a newick, nexus, or nexml
    file.

    There are no constraints on the input trees excepted that they must
    be rooted to perform the analysis. If they are repeated leaves
    (polymorphism), they are automatically removed (several methods
    are implemented, see below).

    Several options are available for analysing the results: the user can
    compute a consensus when several cladograms are optimal, a specific
    character-state testing procedure can be used to test whether each
    character state (i.e., each informative node of a tree) is a synapomorphy
    or an homoplasy, and finally a retention index can be computed to
    obtained the proportion of phylogenetic relationships of characters
    that have been retained in the optimal cladogram.


    Usage:

        agatta analysis <file> [-s -v --analysis=<type> --chartest
                                      --consensus=<type> --parallel=<int>
                                      --pdf=<path> --prefix=<file>
                                      --repetitions=<type> --replicates=<int>
                                      --ri --rosetta=<file> --softpath=<path>
                                      --software=<type> --taxarep1=<path>
                                      --weighting=<type>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --analysis=<type>  Type of tree search analysis between an exact
            branch and bound ('bandb') or an heuristic tree search
            ('heuristic'). The heuristic search is only available through
            PAUP* and TNT, thus it is mandatory to add the flag --software=tnt
            or --software=paup (and the flag --softpath accordingly).
            By default the analysis is in branch and bound below 15 terminals
            and heuristic otherwise.

            --chartest  Test and locates all character states on the cladogram.
            Each character state can be a synapomorphy if the hypothesis is
            accepted or an homoplasy if the hypothesis of sameness is rejected.
            Two output files are writen: prefix.chartest gives all locations
            of the states on the cladogram and if the test is passed or not,
            sorted by state, and prefix.chartest_node gives the same
            information sorted by node.

            --consensus=<type>  Compute a consensus which can be a strict
            consensus ('strict') or a reduced cladistic consensus ('rcc',
            Wilkinson 1994) which is able to detect more common information
            in subtrees. By default, the flag without argument produces a
            strict consensus. The output file is prefix.constrict or
            prefix.rcc.

            --parallel=<type>  Option for choosing if the analysis is made
            using multiprocessing or not. This argument can be:
              - 'not' if the user does not wan to use multiprocessing.
              - 'auto' for automatic detection of the number of cpu.
              - any integer corresponding to the number of cpu allowed.
            By default, the analysis is made in parallel using all available
            cpu. This option can be used in the case of very large character
            tree that can saturate the RAM if too many parallel processing are
            active.

            --pdf=<path>  Compute pdf files to visualise character states on
            the cladogram if --chartest is used. One pdf for each state is
            writen named prefix.character_name.character_state_number.pdf.
            The flag's argument is the path where the user wants to save the
            pdfs.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --repetitions=<type>  The removal of repeated leaves in character
            trees is made using the method of free-paralogy subtree.
            Two algorithms are implemented, the original one from Nelson and
            Ladiges (1996) ('Nelson') for dealing with paralogy in cladistic
            biogeography, and the algorithm of Rineau et al. (2021) ('Rineau')
            designed for all cases of repetitions (not only paralogy).
            If the flag is not used and if repetitions are detected, they are
            automatically removed using Rineau et al. algorithm's.
            The repetition-free character trees are writen in prefix.poly and
            each new tree receives an id, e.g. 1.2 corresponds to the second
            repetition-free subtree computed from the 1st original character
            tree.

            --replicates=<int>  State the number of replicates in case
            of an heuristic tree search. By default the heuristic search is
            launched with 1000 replicates.

            --ri  Compute the retention index of the resulting cladogram and
            a retention index for each character which states the percentage
            of phylogenetic information retained in the optimal cladogram.
            The results are writen in prefix.ri.

            --rosetta=<file>  If the input tree leaves are parts, cladistics
            requires a standardisation step with replacement of parts to
            wholes. It is especially important in cladistic biogeography where
            terminal taxa are replaced by biogeographic areas. The rosetta flag
            replaces leaves of input trees according to a csv file with its
            path given as argument. The csv file is a table with two columns,
            one with the name of the tree leaves and the second with their
            corresponding names to be switched. The results of the
            standardisation are writen in the file prefix.std.

            --softpath=<path>  Path of the software declared in --software.

            --software=<type>  Choose how to perform the three-item analysis.
            The analysis can be performed using the built-in branch and bound
            in Agatta ('agatta'). 'paup' and 'tnt' can be used for branch and
            bound or heuristic search. By default the analysis is made using
            built-in branch and bound. However it works only with very few
            terminals. User should consider to switch the software is the
            analysis time appears to be too long. A prefix.nex file is
            generated if 'paup', and a prefix.tnt file for 'tnt'.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --weighting=<type>  Weighting scheme to use on triplets. The
            type of weighting will change the results of the analysis.
            The following schemes can be used:
                - FW: Fractional weighting from Rineau et al. (2021),
                - FWNL: Fractional weighting from Nelson and Ladiges (1992),
                - MW: Minimal weighting from Wilkinson et al. (2004),
                - AW: Additive weighting : the weight of a triplet in additive
                  weighting corresponds to the number of trees in which the
                  triplet is present,
                - NW: No weighting (all triplets have a weight of 1).
            By default 'FW' is used.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "tripdec":
        print("\033[1;34;48m    tripdec")
        print("""
    Decomposes rooted tree(s) into minimal cladistic statements (triplets
    or three-item statements) stating that two leaves are closer between them
    than to a third. During decomposition the weight of each triplet is
    computed according to a specific weighting scheme. In three-item analysis,
    the weighted triplets are then analysed to compute the cladogram that is
    in agreement with the maximum amout of them.


    Usage:

        agatta tripdec <file> [-s -v --parallel=<int> --prefix=<file>
                                     --taxarep1=<path> --weighting=<type>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --parallel=<type>  Option for choosing if the analysis is made
            using multiprocessing or not. This argument can be:
              - 'not' if the user does not wan to use multiprocessing.
              - 'auto' for automatic detection of the number of cpu.
              - any integer corresponding to the number of cpu allowed.
            By default, the analysis is made in parallel using all available
            cpu. This option can be used in the case of very large character
            tree that can saturate the RAM if too many parallel processing are
            active.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --weighting=<type>  Weighting scheme to use on triplets. The
            type of weighting will change the results of the analysis.
            The following schemes can be used:
                - FW: Fractional weighting from Rineau et al. (2021),
                - FWNL: Fractional weighting from Nelson and Ladiges (1992),
                - MW: Minimal weighting from Wilkinson et al. (2004),
                - AW: Additive weighting : the weight of a triplet in additive
                  weighting corresponds to the number of trees in which the
                  triplet is present,
                - NW: No weighting (all triplets have a weight of 1).
            By default 'FW' is used.


    Output:

        Two files are writen after decomposition. The main one is
        prefix.triplet, a file containing all triplets deduced from the input
        character trees. Each row corresponds to one triplet with its
        weight as an integer (if UW, AW, or NW), or as a fraction and as a
        float otherwise. The second file prefix.taxabloc is a table file with
        the correspondance between leaf identifiers and names given in the
        input trees.
              """)

    elif command == "support":
        print("\033[1;34;48m    support")
        print("""
    The support command gathers several triplet metrics and indices that are
    intended to compare trees between them. The retention index measures the
    amount of cladistic relationships from the input characters retained in
    the optimal cladogram (or a consensus). The triplet distance measures
    the distance in terms of triplets between two trees. The ITRI compares
    in terms of triplets a tree relatively to a reference tree (used to
    measure efficiency of methods using simulations).
    All these metrics handle weighting schemes.

    Usage:

        agatta support <file> <file> [-s -v --index=<type> --prefix=<file>
                                       --taxarep1=<path> --taxarep2=<path>
                                       --weighting=<type>]

        Mandatory parameters:

            The requested files depend of the --index the user wants to
            compute:

                --index=<type>

            Two <file> arguments are requested which represents the path of
            tree files (more on accepted formats here INSERT_URL_HERE)


              Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

            <file>

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --index  prout

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --taxarep2=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --weighting=<type>  Weighting scheme to use on triplets. The
            type of weighting will change the results of the analysis.
            The following schemes can be used:
                - FW: Fractional weighting from Rineau et al. (2021),
                - FWNL: Fractional weighting from Nelson and Ladiges (1992),
                - MW: Minimal weighting from Wilkinson et al. (2004),
                - AW: Additive weighting : the weight of a triplet in additive
                  weighting corresponds to the number of trees in which the
                  triplet is present,
                - NW: No weighting (all triplets have a weight of 1).
            By default 'FW' is used.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "chartest":
        print("\033[1;34;48m    chartest")
        print("""
    Text

    Usage:

        agatta chartest <file> <file> [-s -v --pdf=<path> --prefix=<file>
                                       --taxarep1=<path> --taxarep2=<path>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

            <file>

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --pdf=<path>  Compute pdf files to visualise character states on
            the cladogram if --chartest is used. One pdf for each state is
            writen named prefix.character_name.character_state_number.pdf.
            The flag's argument is the path where the user wants to save the
            pdfs.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --taxarep2=<path>


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "convert":
        print("\033[1;34;48m    convert")
        print("""
    Text

    Usage:

        agatta convert <file> [-s -v --analysis=<type> --filetype=<type> --log
                                      --multiplier=<int> --replicates=<int>
                                      --parallel=<int> --prefix=<file>
                                      --software=<type> --taxarep1=<path>
                                      --weighting=<type>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --analysis=<type>  Type of tree search analysis between an exact
            branch and bound ('bandb') or an heuristic tree search
            ('heuristic'). The heuristic search is only available through
            PAUP* and TNT, thus it is mandatory to add the flag --software=tnt
            or --software=paup (and the flag --softpath accordingly).
            By default the analysis is in branch and bound below 15 terminals
            and heuristic otherwise.

            --filetype=<type>

            --log

            --multiplier=<int>

            --parallel=<type>  Option for choosing if the analysis is made
            using multiprocessing or not. This argument can be:
              - 'not' if the user does not wan to use multiprocessing.
              - 'auto' for automatic detection of the number of cpu.
              - any integer corresponding to the number of cpu allowed.
            By default, the analysis is made in parallel using all available
            cpu. This option can be used in the case of very large character
            tree that can saturate the RAM if too many parallel processing are
            active.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --replicates=<int>  State the number of replicates in case
            of an heuristic tree search. By default the heuristic search is
            launched with 1000 replicates.

            --software=<type>  Choose how to perform the three-item analysis.
            The analysis can be performed using the built-in branch and bound
            in Agatta ('agatta'). 'paup' and 'tnt' can be used for branch and
            bound or heuristic search. By default the analysis is made using
            built-in branch and bound. However it works only with very few
            terminals. User should consider to switch the software is the
            analysis time appears to be too long. A prefix.nex file is
            generated if 'paup', and a prefix.tnt file for 'tnt'.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.

            --weighting=<type>  Weighting scheme to use on triplets. The
            type of weighting will change the results of the analysis.
            The following schemes can be used:
                - FW: Fractional weighting from Rineau et al. (2021),
                - FWNL: Fractional weighting from Nelson and Ladiges (1992),
                - MW: Minimal weighting from Wilkinson et al. (2004),
                - AW: Additive weighting : the weight of a triplet in additive
                  weighting corresponds to the number of trees in which the
                  triplet is present,
                - NW: No weighting (all triplets have a weight of 1).
            By default 'FW' is used.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "fp":
        print("\033[1;34;48m    fp")
        print("""
    Text

    Usage:

        agatta fp <file> [-s -v --prefix=<file> --repetitions=<type>
                                 --taxarep1=<path>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --repetitions=<type>  The removal of repeated leaves in character
            trees is made using the method of free-paralogy subtree.
            Two algorithms are implemented, the original one from Nelson and
            Ladiges (1996) ('Nelson') for dealing with paralogy in cladistic
            biogeography, and the algorithm of Rineau et al. (2021) ('Rineau')
            designed for all cases of repetitions (not only paralogy).
            If the flag is not used and if repetitions are detected, they are
            automatically removed using Rineau et al. algorithm's.
            The repetition-free character trees are writen in prefix.poly and
            each new tree receives an id, e.g. 1.2 corresponds to the second
            repetition-free subtree computed from the 1st original character
            tree.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "consensus":
        print("\033[1;34;48m    consensus")
        print("""
    Text

    Usage:

        agatta consensus <file> [-s -v --consensus=<type> --prefix=<file>
                                       --taxarep1=<path>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --consensus=<type>  Compute a consensus which can be a strict
            consensus ('strict') or a reduced cladistic consensus ('rcc',
            Wilkinson 1994) which is able to detect more common information
            in subtrees. By default, the flag without argument produces a
            strict consensus. The output file is prefix.constrict or
            prefix.rcc.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --taxarep1=<path>  If the user wants to replace identifiers by real
            leaf names in the result files, this flag can be used with a path
            to a csv file with two columns, the first with the identifiers in
            the actual newick strings and the other with the names the user
            wants.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.


    References
              """)

    elif command == "describetree":
        print("\033[1;34;48m    describetree")
        print("""
    Text

    Usage:

        agatta describetree <file> [-s -v --prefix=<file> --showtaxanames]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.

            --showtaxanames  blabla


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "standardisation":
        print("\033[1;34;48m    standardisation")
        print("""
    Text

    Usage:

        agatta standardisation <file> <file> [-s -v --prefix=<file>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

            <file> Path

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    elif command == "hmatrix":
        print("\033[1;34;48m    hmatrix")
        print("""
    Text

    Usage:

        agatta analysis <file> [-s -v --chardec --prefix=<file>]

        Mandatory parameters:

            <file> Path of the file containing the character trees.
            The trees can be encoded for Agatta in a file in several ways:

              - a hierarchical matrix (more on the hierarchical matrix
                format here: INSERT_URL_HERE)

              - A newick file with a single newick tree on each line.
                (see explanations here:
                https://en.wikipedia.org/wiki/Newick_format)

                Example:

                (a,b,(c,d,(e,f)));
                ((c,e),(a,(b,(d,f))));

              - A nexus file (extension in .nex).

              - A nexml file (extension in .nexml).

        Optionnal parameters:

            -s  Silent mode.

            -v  Verbose mode.

            --chardec

            --prefix=<file>  Prefix of all saving files. The complete path can
            be used. By default, the prefix is 'agatta_out' and all files are
            saved in the directory of the first <file>.


    Output:

        Four ouput files are writen all the time when using the analysis
        command in addition to optionnal output files.
        The files are:

          - prefix.log is a log file with all the parameters of the analysis.

          - prefix.triplet is a file with all triplets deduced from the input
            character trees. Each row corresponds to one triplet with its
            weight as a fraction and as a float.

          - prefix.taxabloc is a table file with the correspondance between
            leaf identifiers and names given in the input.

          - prefix.tre is a newick file recording all the optimal cladograms
            found during the analysis.
              """)

    else:
        sys.exit(print("""This command does not exist. Available commands:
              analysis, tripdec, ri, chartest, convert, fp, consensus,
              describetree, standardisation, hmatrix."""))