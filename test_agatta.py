# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@gmail.com

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

from ete3 import Tree
from tqdm import tqdm
from fractions import Fraction
from agatta.search import triplet_check
from agatta.analysis import triplet
from agatta.analysis import standard_tripdec
from agatta.analysis import parallel_tripdec

class Test_analysis:

    def setup_method(self):
        self.nb_try = 20
        self.tree1 = "(a,b,(c,d,(e,f)));"
        self.tree2 = "((c,e),(a,(b,(d,f))));"
        self.character_dict1 = {Tree(self.tree1):1}
        self.character_dict2 = {Tree(self.tree1):1,
                                Tree(self.tree2):2}

        # (a,b,(c,d,(e,f))); in FW
        self.triplet_dict1 = {triplet({"f","e"},{"a"}): Fraction(1, 1),
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
        self.triplet_dict2 = {triplet({"f","e"},{"a"}): Fraction(3, 2),
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


    def test_standard_tripdec_FW(self):

        triplet_dict = standard_tripdec(self.character_dict1, weighting='FW')

        print("test_standard_tripdec_FW")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

        assert triplet_dict == self.triplet_dict1

    def test_parallel_tripdec_FW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='FW')

        print("test_parallel_tripdec_FW")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

        assert triplet_dict == self.triplet_dict1

    def test_standard_tripdec_FWNL(self):

        triplet_dict = standard_tripdec(self.character_dict1, weighting='FWNL')

        print("test_standard_tripdec_FWNL")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

        assert triplet_dict == self.triplet_dict2

    def test_parallel_tripdec_FWNL(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='FWNL')

        print("test_parallel_tripdec_FWNL")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

        assert triplet_dict == self.triplet_dict2

    def test_standard_tripdec_NW(self):

        triplet_dict = standard_tripdec(self.character_dict1, weighting='NW')

        print("test_standard_tripdec_NW")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

        assert triplet_dict == self.triplet_dictNW

    def test_parallel_tripdec_NW(self):

        triplet_dict = parallel_tripdec(self.character_dict1, weighting='NW')

        print("test_parallel_tripdec_NW")
        for trip, weight in triplet_dict.items():
            print(str(trip) + ":    " + str(weight))

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

                    for t, w in triplet_dict_s.items():
                        if not t in triplet_dict_p:
                            print(t)
                            print("standard tripdec")
                        if not w == 1:
                            print(t)
                            print("standard tripdec weights")

                    for t, w in triplet_dict_p.items():
                        if not t in triplet_dict_s:
                            print(t)
                            print("parallel tripdec")
                        if not w == 1:
                            print(t)
                            print(rand_tree.write(format=9))
                            print(rand_tree)
                            print("parallel tripdec weights")
                            break

                    print(len(triplet_dict_s))
                    print(len(triplet_dict_p))
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

                    for t, w in triplet_dict_s.items():
                        if not triplet_check(t, cladogram):
                            print(t)
                            print("standard tripdec")

                    for t, w in triplet_dict_p.items():
                        if triplet_check(t, cladogram):
                            print(t)
                            print("parallel tripdec")

                    success = False
                    break

        assert success == True
