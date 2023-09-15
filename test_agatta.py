# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from ete3 import Tree
from tqdm import tqdm
from fractions import Fraction
from agatta.search import triplet_check
from agatta.ini import character_extraction
from agatta.ini import hmatrix
from agatta.ini import standardisation
from agatta.interpret import constrict
from agatta.interpret import rcc
from agatta.interpret import RI
from agatta.interpret import triplet_distance
from agatta.interpret import character_states_test
from agatta.interpret import describe_forest
from agatta.interpret import NRI
from agatta.out import triplet_nexus_file
from agatta.out import triplet_tnt_file
from agatta.analysis import triplet
from agatta.analysis import standard_tripdec
from agatta.analysis import parallel_tripdec
from agatta.analysis import del_replications

class Test_analysis:

    def setup_method(self):
        self.nb_try = 1
        self.tree1 = "(a,b,(c,d,(e,f)));"
        self.tree2 = "((c,e),(a,(b,(d,f))));"
        self.tree3 = "(a,b,(c,e,(d,f)));"
        self.tree4 = "(a,b,(c,d,e,f));"
        self.pectree = ("(a,(b,(c,(d,(e,(f,(g,(h,(i,(j,(k,(l,(m,("
                        + "n,(o,(p,(q,(r,(s,(t,(u,(v,(w,(x,(y,z)))))))"
                        + "))))))))))))))))));")
        self.treerep = "(a,(b,(c,(d,(e,(f,a,(g,(h,i))),(j,(k,(l,i))))))));"
        self.rcc1 = "(a,(b,(c,(d,(e,(f,(g,h)))))));"
        self.rcc2 = "((a,(b,f)),(c,((d,g),(h,e))));"
        self.character_dict1 = {Tree(self.tree1):1}
        self.character_dict2 = {Tree(self.tree1):1,
                                Tree(self.tree2):2}
        self.treelist = [Tree(self.tree1), Tree(self.tree3)]
        self.pectinate_tree = {Tree(self.pectree): 1}
        self.std = "(France,France,(Italy,(Greece,Spain)));"
        self.hmatrix = ("(Diceras,(Hippurites,(Radiolites,Titanosarcolites),"
                       + "(Titanosarcolites,Clinocaprina,Vaccinites)));")
        self.treerep = "(a,(b,(c,(d,(e,(f,a,(g,(h,i))),(j,(k,(l,i))))))));"
        self.noreptree = ['(b,c,d,e,j,k,l,(f,a,(g,(h,i))));',
                          '(a,b,c,d,e,f,g,h,(j,(k,(l,i))));',
                          '(b,(c,(d,(e,f,a,g,h,i,j,k,l))));']
        self.noreptreeFPS = ['(a,(b,(c,(d,(j,(k,(l,i)))))));',
                             '(b,(c,(d,(f,a,(g,(h,i))))));']

        self.syn_dict2 = {'a': {'accepted': [], 'rejected': []},
         'b': {'accepted': [], 'rejected': []},
         'c': {'accepted': [], 'rejected': []},
         'e': {'accepted': [], 'rejected': []},
         'd': {'accepted': [], 'rejected': []},
         'f': {'accepted': [], 'rejected': []},
         '2': {'accepted': [], 'rejected': []},
         '1': {'accepted': [['#1.1']], 'rejected': []},
         'root': {'accepted': [], 'rejected': []}}


        # (a,b,(c,d,(e,f))); in FW
        self.triplet_dictFW = {triplet({"f","e"},{"a"}): Fraction(1, 1),
                              triplet({"f","e"},{"b"}): Fraction(1, 1),
                              triplet({"f","e"},{"c"}): Fraction(1, 1),
                              triplet({"f","e"},{"d"}): Fraction(1, 1),
                              triplet({"c","d"},{"a"}): Fraction(3, 5),
                              triplet({"c","e"},{"a"}): Fraction(3, 5),
                              triplet({"c","f"},{"a"}): Fraction(3, 5),
                              triplet({"d","e"},{"a"}): Fraction(3, 5),
                              triplet({"d","f"},{"a"}): Fraction(3, 5),
                              triplet({"c","d"},{"b"}): Fraction(3, 5),
                              triplet({"c","e"},{"b"}): Fraction(3, 5),
                              triplet({"c","f"},{"b"}): Fraction(3, 5),
                              triplet({"d","e"},{"b"}): Fraction(3, 5),
                              triplet({"d","f"},{"b"}): Fraction(3, 5)}

        # (a,b,(c,d,(e,f))); in FWNL
        self.triplet_dictFWNL = {triplet({"f","e"},{"a"}): Fraction(3, 2),
                                triplet({"f","e"},{"b"}): Fraction(3, 2),
                                triplet({"f","e"},{"c"}): Fraction(1, 1),
                                triplet({"f","e"},{"d"}): Fraction(1, 1),
                                triplet({"c","d"},{"a"}): Fraction(1, 2),
                                triplet({"c","e"},{"a"}): Fraction(1, 2),
                                triplet({"c","f"},{"a"}): Fraction(1, 2),
                                triplet({"d","e"},{"a"}): Fraction(1, 2),
                                triplet({"d","f"},{"a"}): Fraction(1, 2),
                                triplet({"c","d"},{"b"}): Fraction(1, 2),
                                triplet({"c","e"},{"b"}): Fraction(1, 2),
                                triplet({"c","f"},{"b"}): Fraction(1, 2),
                                triplet({"d","e"},{"b"}): Fraction(1, 2),
                                triplet({"d","f"},{"b"}): Fraction(1, 2)}

        # (a,b,(c,d,(e,f))); in MW
        self.triplet_dictMW = {triplet({"f","e"},{"a"}): Fraction(1, 5),
                                triplet({"f","e"},{"b"}): Fraction(1, 5),
                                triplet({"f","e"},{"c"}): Fraction(1, 5),
                                triplet({"f","e"},{"d"}): Fraction(1, 5),
                                triplet({"c","d"},{"a"}): Fraction(1, 5),
                                triplet({"c","e"},{"a"}): Fraction(1, 5),
                                triplet({"c","f"},{"a"}): Fraction(1, 5),
                                triplet({"d","e"},{"a"}): Fraction(1, 5),
                                triplet({"d","f"},{"a"}): Fraction(1, 5),
                                triplet({"c","d"},{"b"}): Fraction(1, 5),
                                triplet({"c","e"},{"b"}): Fraction(1, 5),
                                triplet({"c","f"},{"b"}): Fraction(1, 5),
                                triplet({"d","e"},{"b"}): Fraction(1, 5),
                                triplet({"d","f"},{"b"}): Fraction(1, 5)}

        # (a,b,(c,d,(e,f))); in AW
        self.triplet_dictAW = {triplet({"f","e"},{"a"}): 1,
                              triplet({"f","e"},{"b"}): 1,
                              triplet({"f","e"},{"c"}): 1,
                              triplet({"f","e"},{"d"}): 1,
                              triplet({"c","d"},{"a"}): 1,
                              triplet({"c","e"},{"a"}): 1,
                              triplet({"c","f"},{"a"}): 1,
                              triplet({"d","e"},{"a"}): 1,
                              triplet({"d","f"},{"a"}): 1,
                              triplet({"c","d"},{"b"}): 1,
                              triplet({"c","e"},{"b"}): 1,
                              triplet({"c","f"},{"b"}): 1,
                              triplet({"d","e"},{"b"}): 1,
                              triplet({"d","f"},{"b"}): 1}

        # (a,b,(c,d,(e,f))); in UW
        self.triplet_dictUW = {triplet({"f","e"},{"a"}): 2,
                              triplet({"f","e"},{"b"}): 2,
                              triplet({"f","e"},{"c"}): 1,
                              triplet({"f","e"},{"d"}): 1,
                              triplet({"c","d"},{"a"}): 1,
                              triplet({"c","e"},{"a"}): 1,
                              triplet({"c","f"},{"a"}): 1,
                              triplet({"d","e"},{"a"}): 1,
                              triplet({"d","f"},{"a"}): 1,
                              triplet({"c","d"},{"b"}): 1,
                              triplet({"c","e"},{"b"}): 1,
                              triplet({"c","f"},{"b"}): 1,
                              triplet({"d","e"},{"b"}): 1,
                              triplet({"d","f"},{"b"}): 1}

        # (a,b,(c,d,(e,f))); in NW
        self.triplet_dictNW = {triplet({"f","e"},{"a"}): 1,
                              triplet({"f","e"},{"b"}): 1,
                              triplet({"f","e"},{"c"}): 1,
                              triplet({"f","e"},{"d"}): 1,
                              triplet({"c","d"},{"a"}): 1,
                              triplet({"c","e"},{"a"}): 1,
                              triplet({"c","f"},{"a"}): 1,
                              triplet({"d","e"},{"a"}): 1,
                              triplet({"d","f"},{"a"}): 1,
                              triplet({"c","d"},{"b"}): 1,
                              triplet({"c","e"},{"b"}): 1,
                              triplet({"c","f"},{"b"}): 1,
                              triplet({"d","e"},{"b"}): 1,
                              triplet({"d","f"},{"b"}): 1}

        # (a,b,(c,d,(e,f))); for bandb test
        self.triplet_dict_outfile = {triplet({"f","e"},{"a"}): 1,
                                      triplet({"f","e"},{"b"}): 1,
                                      triplet({"f","e"},{"c"}): 1,
                                      triplet({"f","e"},{"d"}): 1,
                                      triplet({"c","d"},{"a"}): 1,
                                      triplet({"c","e"},{"a"}): 1,
                                      triplet({"c","f"},{"a"}): 1,
                                      triplet({"d","e"},{"a"}): 2,
                                      triplet({"d","f"},{"a"}): 1,
                                      triplet({"c","d"},{"b"}): 1,
                                      triplet({"c","e"},{"b"}): 1,
                                      triplet({"c","f"},{"b"}): 1,
                                      triplet({"d","e"},{"b"}): 1,
                                      triplet({"d","f"},{"b"}): 1}

        self.pie_percentages1 = {'1': [100], '2': [100], '3': [100]}

    def test_character_extraction(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        infile = str(folder.join("input.txt"))

        with open(infile, "w") as treefile:
            treefile.write(self.tree1 + "\n")

        character_dict = character_extraction(infile)
        result_tree = list(character_dict.keys())[0]
        result_tree.ladderize()
        result_tree = result_tree.write(format=9)

        assert result_tree == self.tree1


    def test_standard_tripdec_FW(self):

        triplet_dict = standard_tripdec(self.character_dict1, weighting='FW')

        assert triplet_dict == self.triplet_dictFW


    def test_parallel_tripdec_FW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='FW')

        assert triplet_dict == self.triplet_dictFW


    def test_standard_tripdec_FWNL(self):

        triplet_dict = standard_tripdec(self.character_dict1, weighting='FWNL')

        assert triplet_dict == self.triplet_dictFWNL


    def test_parallel_tripdec_FWNL(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='FWNL')

        assert triplet_dict == self.triplet_dictFWNL


    def test_parallel_tripdec_MW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='MW')

        assert triplet_dict == self.triplet_dictMW


    def test_parallel_tripdec_UW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='UW')

        assert triplet_dict == self.triplet_dictUW


    def test_parallel_tripdec_AW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='AW')

        assert triplet_dict == self.triplet_dictAW


    def test_parallel_tripdec_NW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='NW')

        assert triplet_dict == self.triplet_dictNW


    def test_identical_parallel_standard_FW(self):

        triplet_dict_s = standard_tripdec(self.character_dict2, weighting='FW')
        triplet_dict_p = parallel_tripdec(self.character_dict2, weighting='FW')

        assert triplet_dict_s == triplet_dict_p


    def test_identical_parallel_standard_random(self):

        for weight in ["FW","FWNL","UW","MW","AW","NW"]:
            character_dict = dict()
            success = True

            for _ in tqdm(range(self.nb_try), desc=weight):
                for i in range(2):
                    rand_tree = Tree()
                    rand_tree.populate(15)
                    character_dict[rand_tree] = 0

                triplet_dict_s = dict(standard_tripdec(character_dict,
                                                  weighting=weight))
                triplet_dict_p = dict(parallel_tripdec(character_dict,
                                                  weighting=weight))

                if not triplet_dict_s == triplet_dict_p:
                    success = False
                    break

        assert success == True


    def test_triplet_in_tree(self):

        success = True

        for weight in ["FW","FWNL","UW","MW","AW","NW"]:
            for _ in tqdm(range(self.nb_try), desc=weight):
                character_dict = dict()
                cladogram = Tree()
                cladogram.populate(15)
                character_dict[cladogram] = 0

                triplet_dict_s = dict(standard_tripdec(character_dict,
                                                  weighting=weight))
                triplet_dict_p = dict(parallel_tripdec(character_dict,
                                                  weighting=weight))

                if not triplet_dict_s == triplet_dict_p:
                    success = False

        assert success == True


    def test_pectinate_tree_FW(self):
        """
        check if fully pectinate tree always produces triplet weights of 1
        with FW.
        """

        triplet_dict_s = dict(standard_tripdec(self.pectinate_tree,
                                          weighting="FW"))
        triplet_dict_p = dict(parallel_tripdec(self.pectinate_tree,
                                          weighting="FW"))

        assert (all(value == 1 for value in triplet_dict_s.values()) and
                all(value == 1 for value in triplet_dict_p.values()))

    def test_standardisation(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        infile = str(folder.join("std_input.tre"))
        biogeo_tab = str(folder.join("biogeo_tab.csv"))
        outfile = str(folder.join("std_output.tre"))

        with open(infile, "w") as treefile:
            treefile.write(self.tree1 + "\n")

        with open(biogeo_tab, "w") as treefile:
            treefile.write("a;France\n")
            treefile.write("b;France\n")
            treefile.write("c;Czechia\n")
            treefile.write("c;Germany\n")
            treefile.write("d;Italy\n")
            treefile.write("e;Greece\n")
            treefile.write("f;Spain\n")

        areagram_dict = standardisation(infile,
                                        biogeo_tab,
                                        outfile,
                                        verbose=False)

        result_tree = list(areagram_dict.keys())[0]
        result_tree.ladderize()
        result_tree = result_tree.write(format=9)

        assert result_tree == self.std


    def test_hmatrix(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        infile = str(folder.join("input.txt"))

        with open(infile, "w") as treefile:
            treefile.write(";(0,(1,(2),(3)))\n")
            treefile.write("Diceras;0\n")
            treefile.write("Hippurites;1\n")
            treefile.write("Radiolites;2\n")
            treefile.write("Titanosarcolites;2,3\n")
            treefile.write("Clinocaprina;3\n")
            treefile.write("Vaccinites;3\n")

        character_dict = hmatrix(infile)

        result_tree = list(character_dict.keys())[0]
        result_tree.ladderize()
        result_tree = result_tree.write(format=9)

        assert result_tree == self.hmatrix


    def test_constrict(self):

        constree = constrict(self.treelist)

        assert constree.write(format=9) == self.tree4


    def test_rcc(self):
        """
        Test example of figure 2 p350 in Wilkinson 1994

        """

        profile = rcc([Tree(self.rcc1), Tree(self.rcc2)])

        assert len(profile) == 5


    def test_retention_index(self):

        RI_char_dict = RI(self.character_dict1, self.character_dict2, 
                          verbose=True)

        assert RI_char_dict == {'1': Fraction(1, 1), '2': Fraction(1, 5),
                      'Total': Fraction(7, 15)}


    def test_triplet_distance(self):

        Precision, Recall, Fscore = triplet_distance(Tree(self.tree1),
                                                     Tree(self.tree2),
                                                     prefix=False)
        
        assert [Precision, Recall, Fscore] == [0.2, 0.24, 0.2181818181818182]


    def test_charstate(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        prefix = str(folder.join("input"))

        results_dict, syn_dict = character_states_test({Tree(self.tree3):1},
                                                       {Tree(self.tree4):1},
                                                       prefix)

        assert syn_dict == self.syn_dict2


    def test_describe_forest(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        prefix = str(folder.join("input"))

        describe_forest(self.character_dict1, prefix)

        with open(prefix + ".dt", "r") as describefile:
            for line in describefile:
                if line.startswith("Number of nodes"):

                    assert line.strip().endswith("9")


    def test_triplet_nexus_file(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        prefix = str(folder.join("input"))


        triplet_nexus_file(self.triplet_dict_outfile,
                            self.character_dict1,
                            weighting="FW",
                            prefix=prefix,
                            analysis="bandb")

        wts = False

        with open(prefix + ".nex", "r") as describefile:
            for line in describefile:
                if line.startswith("a "):
                    a_line = line.strip().split(" ")[1]
                if line.startswith("d "):
                    d_line = line.strip().split(" ")[1]
                if line.startswith("e "):
                    e_line = line.strip().split(" ")[1]
                if line.startswith("wts"):
                    wts = True
                if wts and line[1] == ":":
                    w_line = [w[0] for w in line.strip().split(" ")]
                    w_number = w_line.index('2')
                    wts = False

        assert (a_line[w_number] == '0'
                and d_line[w_number] == '1'
                and e_line[w_number] == '1')


    def test_triplet_tnt_file(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        prefix = str(folder.join("input"))


        triplet_tnt_file(self.triplet_dict_outfile,
                            self.character_dict1,
                            weighting="FW",
                            prefix=prefix,
                            analysis="bandb")

        with open(prefix + ".tnt", "r") as describefile:
            for line in describefile:
                if line.startswith("a "):
                    a_line = line.strip().split(" ")[1]
                if line.startswith("d "):
                    d_line = line.strip().split(" ")[1]
                if line.startswith("e "):
                    e_line = line.strip().split(" ")[1]
                if line.startswith("ccode"):
                    w_line = [int(w[1:]) for w in line.strip().split(" ")[1:]
                              if w.startswith('/')]
                    w_number = w_line.index(1000)

        assert (a_line[w_number] == '0'
                and d_line[w_number] == '1'
                and e_line[w_number] == '1')


    # def test_bandb(self):

    #     optimal_score, optimal_tree_list = bandb(["a","b","c","d","e","f"],
    #                                              self.triplet_dict_outfile)

    #     constree = constrict(optimal_tree_list)
    #     constree.ladderize()

    #     assert constree.write(format=9) == self.tree1


    def test_del_replications_TMS(self):

        subtrees = del_replications(Tree(self.treerep))

        assert [t.write(format=9) for t in subtrees] == self.noreptree


    def test_del_replications_FPS(self):

        subtrees = del_replications(Tree(self.treerep), method="FPS")

        assert [t.write(format=9) for t in subtrees] == self.noreptreeFPS


    def test_nri(self, tmpdir):

        folder = tmpdir.mkdir("subdir")
        infile = str(folder.join("input.hmatrix"))
        cladogramfile = str(folder.join("cladogram.txt"))


        with open(infile, "w") as treefile:
            treefile.write(";(0,(1,(2),(3)))\n")
            treefile.write("Diceras;0\n")
            treefile.write("Hippurites;1\n")
            treefile.write("Radiolites;2\n")
            treefile.write("Titanosarcolites;2,3\n")
            treefile.write("Clinocaprina;3\n")
            treefile.write("Vaccinites;3\n")

        with open(cladogramfile, "w") as treefile2:
            treefile2.write("(Diceras,(Hippurites,(Radiolites,Vaccinites),(Titanosarcolites,Clinocaprina)));\n")

        pie_percentages1 = NRI(cladogramfile, infile, 
            taxarep1=False, 
            taxarep2=False, 
            prefix=str(folder.join("rnri")), 
            rnri_codes=False, 
            weighting='FW', 
            polymethod='TMS', 
            totaltree=True, 
            bubble_size=0.05, 
            rescaled=True,
            pdf_files=False)

        assert pie_percentages1 == self.pie_percentages1