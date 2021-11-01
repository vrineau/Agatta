# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:42 2020

@author: Valentin Rineau
"""

from fractions import Fraction
from itertools import combinations
#from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from functools import partial
import warnings
import pickle
import uuid
import os

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
    from ete3 import Tree
    

class triplet:
    __slots__ = ("in_taxa", "out_taxa", "FW")
    
    def __init__(self, in_taxa = [], out_taxa = str(), FW = Fraction(), parent_triplets = ()):
        self.in_taxa = in_taxa
        self.out_taxa = out_taxa # taxon out
        self.FW = FW # Liste contenant la fraction de FW, avec d'abord le numérateur puis le dénominateur
        self.parent_triplets = parent_triplets # Liste contenant les triplets parents
            
    def __eq__(self, other_triplet):
        return (set(self.in_taxa) == set(other_triplet.in_taxa)) and (self.out_taxa == other_triplet.out_taxa)

    def __call__(self, in_taxa = [], out_taxa = str(), FW = Fraction(), parent_triplets = ()):
        return self
    
    def __del__(self):
        del self

    def __repr__(self):
        return "({}({},{})) - FW: {}".format(self.out_taxa, self.in_taxa[0], self.in_taxa[1], self.FW)

class id_triplet(triplet):
    def __init__(self, in_taxa = set(), out_taxa = set(), FW = Fraction(), parent_triplets = set(), primary_parent_triplets = set(), gen_level = int()):
        self.in_taxa = set(in_taxa)
        self.out_taxa = set(out_taxa) # taxon out
        self.FW = FW # Liste contenant la fraction de FW, avec d'abord le numérateur puis le dénominateur
        self.parent_triplets = parent_triplets # Liste contenant les triplets parents (parenté directe)
        self.primary_parent_triplets = primary_parent_triplets # Liste contenant les triplets parents primaires
        self.gen_level = gen_level
            
    def __eq__(self, other_triplet):
        return (set(self.in_taxa) == set(other_triplet.in_taxa)) and (self.out_taxa == other_triplet.out_taxa) and set(self.parent_triplets) == set(other_triplet.parent_triplets)

    def __hash__(self):
        return hash(list(self.in_taxa)[0]) ^ hash(list(self.in_taxa)[1]) ^ hash(list(self.out_taxa)[0]) ^ (hash(tuple(self.parent_triplets)))

    def __repr__(self):
        return "({}({},{}))".format(list(self.out_taxa)[0], list(self.in_taxa)[0], list(self.in_taxa)[1])


def tripdecFW(triplet_output=dict(), total_taxaset=False, character=Tree()):
    """
    Essai de calcul des triplets et des poids en FW en récursif

    """
    
    if not total_taxaset:
        total_taxaset = set(character.get_leaf_names())
        
    #génération du dictionnaire de triplets
    tree_tripdic = defaultdict(int) #dictionnaire de triplets de l'arbre
    nodaltripletdict = defaultdict(int) #dictionnaire de triplets de la composante
    children_generator = (child_node for child_node in character.get_children() if child_node.is_leaf() == False) #générateur de tous les noeuds informatifs directement descendants 
    
    #récolte des triplets des noeuds enfants
    for child_node in children_generator: #donc pour chaque couple de noeuds informatifs parent-enfant
        newtripletdict = tripdecFW(defaultdict(int), total_taxaset, character=child_node)
        
        for trip, FW in newtripletdict.items():
            tree_tripdic[trip] += FW

    #calcul des triplets pour ce noeud si non racine
    if character.is_root(): #si noeud racine (l'analyse est terminée)

        for trip in tree_tripdic.keys():
            if trip in triplet_output: #nouveaux triplets
                triplet_output[trip] += tree_tripdic[trip]  
            else: #incrémentation de triplets préexistants
                triplet_output[trip] = tree_tripdic[trip]
    
        return triplet_output #renvoie un dictionnaire de triplets / FW si pas de multiprocessing
        
    else: #sinon
        taxa_in = set(character.get_leaf_names()) #obtenir la liste des taxons d'un noeud - IN
        taxa_out = total_taxaset - taxa_in #obtenir la liste des taxons hors d'un noeud - OUT
        totalnodeweight = len(taxa_out)*(len(taxa_in)-1)
        
        for taxa_in1, taxa_in2 in combinations(taxa_in, 2): #obtenir les combinaisons in_taxa
            for singleout in taxa_out: #construction des triplets et ajout au set du noeud
                if not id_triplet({taxa_in1, taxa_in2}, {singleout}) in tree_tripdic:
                    nodaltripletdict[id_triplet({taxa_in1, taxa_in2}, {singleout})] = 0
                
        for trip in nodaltripletdict.keys():
            tree_tripdic[trip] += Fraction(totalnodeweight, len(nodaltripletdict))

        return tree_tripdic


def tripdec(weighting, character):
    """
        prend un arbre et renvoie une liste de triplets
    """
    #calcul des triplets par caractère
    
    #start = timer()
    #print("construction du dictionnaire de triplets-poids")
    
    triplet_output=dict()
    cardinal_character1 = set(Tree.get_leaf_names(character)) #obtenir la liste des taxons de l'arbre
    internal_nodes1 = (node for node in character.traverse(strategy="preorder") if node.is_leaf() == False and node.is_root() == False)
    tree_tripdic = defaultdict(int)
    card = len(cardinal_character1)
    MW_node = Fraction(6,card*(card-1))
    #calculer les triplets par noeud => construire une liste de triplets
    for node in internal_nodes1:
        taxa_in = set(Tree.get_leaf_names(node)) #obtenir la liste des taxons d'un noeud - IN
        taxa_out = cardinal_character1 - taxa_in #obtenir la liste des taxons hors d'un noeud - OUT    
        FW_node = Fraction(2, len(taxa_in))
        
        for taxa_in1, taxa_in2 in combinations(taxa_in, 2): #obtenir les combinaisons in_taxa
            for singleout in taxa_out: #construction des triplets et ajout au set du noeud                    
                if weighting == "FWNL":
                    tree_tripdic[id_triplet({taxa_in1, taxa_in2}, {singleout})] += FW_node
                elif weighting == "UW":
                    tree_tripdic[id_triplet({taxa_in1, taxa_in2}, {singleout})] += 1
                if weighting == "MW":
                    tree_tripdic[id_triplet({taxa_in1, taxa_in2}, {singleout})] = MW_node
                elif weighting == "AW" or weighting == "NW":
                    tree_tripdic[id_triplet({taxa_in1, taxa_in2}, {singleout})] = 1

    #end = timer()
    #print(f'elapsed time (construction dictionnaire de triplets-poids): {end - start}')
    #start = timer()
    #print("fusion des dictionnaires")

    if weighting == "NW": #si NW
            
        triplet_output.update(tree_tripdic)
            
        return triplet_output
        
    else:
        for trip in tree_tripdic.keys():
            if trip in triplet_output: #nouveaux triplets
                triplet_output[trip] += tree_tripdic[trip]  
            else: #incrémentation de triplets préexistants
                triplet_output[trip] = tree_tripdic[trip]
    
        return triplet_output #renvoie un dictionnaire de triplets / FW si pas de multiprocessing


def tripdecFW_parallel(triplet_output=dict(), total_taxaset=False, pickle_file=False, taxa_dict=dict(), character=Tree()):
    """
    Calcul des triplets et des poids en FW en récursif, en parallèle (multiprocessing)

    """
    #first recursive loop
    if not total_taxaset:
        total_taxaset = set(character.get_leaf_names())
        pickle_name = str(uuid.uuid4()) + "_agatta.pickle"
        pickle_file = open(pickle_name, 'wb')
    
    i = 0
            
    #génération du dictionnaire de triplets
    nodaltripletdict = defaultdict(int) #dictionnaire de triplets de la composante
    children_generator = (child_node for child_node in character.get_children() if child_node.is_leaf() == False) #générateur de tous les noeuds informatifs directement descendants 
    
    #dump des triplets des noeuds enfants et addition du compteur de triplets
    for child_node in children_generator: #donc pour chaque couple de noeuds informatifs parent-enfant
        i += tripdecFW_parallel(defaultdict(int), total_taxaset, pickle_file=pickle_file, character=child_node)

    #si noeud racine, analyse terminée
    if character.is_root():
        
        pickle_file.close()
        
        return [pickle_name, i]
    
    #calcul des triplets pour ce noeud si non racine
    else: #sinon
        taxa_in = set(character.get_leaf_names()) #obtenir la liste des taxons d'un noeud - IN
        taxa_out = total_taxaset - taxa_in #obtenir la liste des taxons hors d'un noeud - OUT
        totalnodeweight = len(taxa_out)*(len(taxa_in)-1)
        
        for taxa_in1, taxa_in2 in combinations(taxa_in, 2): #obtenir les combinaisons in_taxa
            for singleout in taxa_out: #construction des triplets et ajout au set du noeud
                
                #on vérifie si le triplet est original ou s'il existe déjà
                original_triplet = True

                for child_node in (child_node for child_node in character.get_children() if child_node.is_leaf() == False):
                    if taxa_in1 in child_node.get_leaf_names() and taxa_in2 in child_node.get_leaf_names():
                        original_triplet = False
                
                if original_triplet:
                    nodaltripletdict[(taxa_in1, taxa_in2, singleout)] = 0
        
        #dump trip in pickle
        for trip in nodaltripletdict.keys():
            FW = Fraction(totalnodeweight, len(nodaltripletdict))
            pickle.dump((taxa_dict[trip[0]], taxa_dict[trip[1]], taxa_dict[trip[2]], FW), pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
            i += 1

        return i


def tripdec_parallel(weighting, taxa_dict, character):   
    """
        prend un arbre et renvoie une liste de triplets
    """
    #calcul des triplets par caractère
    
    #start = timer()
    #print("construction du dictionnaire de triplets-poids")
    
    pickle_name = str(uuid.uuid4()) + "_agatta.pickle"
    i = 0
    
    cardinal_character1 = set(Tree.get_leaf_names(character)) #obtenir la liste des taxons de l'arbre
    internal_nodes1 = (node for node in character.traverse(strategy="preorder") if node.is_leaf() == False and node.is_root() == False)
    card = len(cardinal_character1)
    MW_node = Fraction(6,card*(card-1))
    #calculer les triplets par noeud => construire une liste de triplets
    with open(pickle_name, 'wb') as pickle_file:
        
        for node in internal_nodes1:
            taxa_in = set(Tree.get_leaf_names(node)) #obtenir la liste des taxons d'un noeud - IN
            taxa_out = cardinal_character1 - taxa_in #obtenir la liste des taxons hors d'un noeud - OUT    
            FW_node = Fraction(2, len(taxa_in))
            
            for taxa_in1, taxa_in2 in combinations(taxa_in, 2): #obtenir les combinaisons in_taxa
                for singleout in taxa_out: #construction des triplets et ajout au set du noeud                    
                    if weighting == "FWNL":
                        pickle.dump((taxa_dict[taxa_in1], taxa_dict[taxa_in2], taxa_dict[singleout], FW_node), pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
                    elif weighting == "MW":
                        pickle.dump((taxa_dict[taxa_in1], taxa_dict[taxa_in2], taxa_dict[singleout], MW_node), pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
                    else: # UW AW ou NW
                        pickle.dump((taxa_dict[taxa_in1], taxa_dict[taxa_in2], taxa_dict[singleout], 1), pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
                        
                    i += 1 #compteur de triplets

    #end = timer()
    #print(f'elapsed time (construction dictionnaire de triplets-poids): {end - start}')
            
    return [pickle_name, i]


def standard_tripdec(character_dict, weighting):
    """
    Prend un tuple d'arbres et compute toutes les listes de triplets en parallèle
    Combine les listes en une seule
    Renvoie la liste
    
    Fonctionne pour tous types de pondérations

    """
    
    #computation des dictionnaires de triplets pour chaque arbre en parallèle
    print("starting {} computations".format(weighting))
    #start = timer()
    
    triplet_output = dict() 
    
    if weighting == "FW":
        for treedec in character_dict.keys():
            triplet_output = tripdecFW(triplet_output, False, character=treedec)
 
    else: #MW, AW, NW
        for treedec in character_dict.keys():
            triplet_output = tripdec(weighting, treedec)
 
    #end = timer()
    #print(f'elapsed time (triplet decomposition and weighting): {end - start}')
    
    return triplet_output


def parallel_tripdec(character_dict, weighting, ncpu="auto"):
    """
    Prend un tuple d'arbres et compute toutes les listes de triplets en parallèle
    Combine les listes en une seule
    Renvoie la liste
    
    Fonctionne pour tous types de pondérations

    """
    
    def taxa_to_numbers(taxa_list): 
        """
            prend un ensemble de taxons et renvoie les taxons avec un identifiant 
            (nombre associé)
        """
        taxa_dict = {}
        taxa_cptr = 0
        for taxa1 in taxa_list: #pour chaque taxon
            taxa_dict[taxa1] = taxa_cptr #integer labels
            taxa_cptr += 1
            
        return taxa_dict
    
    #nombre de coeurs    
    if ncpu == "auto":
        ncpu = cpu_count()
    else:
        ncpu = int(ncpu)

    #computation des dictionnaires de triplets pour chaque arbre en parallèle
    #start = timer()
    print("starting {} computations on {} cores".format(weighting, str(ncpu)))
    
    treetuple = tuple(character_dict) #dictionnaire en tuple pour multiprocessing
    
    #création du dictionnaire 
    taxanames = list(set([x for l in [character.get_leaf_names() for character in character_dict] for x in l]))
    taxa_dict = taxa_to_numbers(taxanames)
    
    if weighting == "FW":
        with Pool(ncpu) as pool: #activation du multiprocessing
            func = partial(tripdecFW_parallel, dict(), False, False, taxa_dict)
            pickle_double = pool.map(func, treetuple) #liste de dictionnaires de triplets

    else:
        with Pool(ncpu) as pool: #activation du multiprocessing
            func = partial(tripdec_parallel, weighting, taxa_dict)
            pickle_double = pool.map(func, treetuple) #liste de dictionnaires de triplets
     
    #end = timer()
    #print(f'elapsed time (triplet decomposition and weighting): {end - start}')
    
    #étape d'extraction des dictionnaires du pickle et combinaison

    #start = timer()
    #print("extraction des dictionnaires du pickle et combinaison")  
        
    shared_dict = defaultdict(int)
    
    if weighting == "FW" or weighting == "FWNL" or weighting == "UW": #poids par composante
    
        for pd in pickle_double:        
            cptr = 0
            with open(pd[0], 'rb') as pickle_file:
                            
                while cptr != pd[1]:
                            
                    triplist = pickle.load(pickle_file)
                    trip = id_triplet({triplist[0], triplist[1]}, {triplist[2]})                       
                    shared_dict[trip] += triplist[3]
                    cptr += 1
                    
    if weighting == "AW" or weighting == "MW": #poids par caractere
    
        for pd in pickle_double:
            cptr = 0
            with open(pd[0], 'rb') as pickle_file: #pour chaque caractere
                
                temp_dict = defaultdict(int)
                while cptr != pd[1]:
                            
                    triplist = pickle.load(pickle_file)
                    trip = id_triplet({triplist[0], triplist[1]}, {triplist[2]})
                    
                    temp_dict[trip] = triplist[3]
                    cptr += 1
                    
                for triptemp, FW in temp_dict.items():
                    shared_dict[triptemp] += FW
                    
    if weighting == "NW": #poids pour tous les caracteres
    
        for pd in pickle_double:
            
            cptr = 0
            
            with open(pd[0], 'rb') as pickle_file:
                            
                while cptr != pd[1]:
                            
                    triplist = pickle.load(pickle_file)
                    trip = id_triplet({triplist[0], triplist[1]}, {triplist[2]})                       
                    shared_dict[trip] = triplist[3]
                    cptr += 1
    
    #supression des pickle
    for pd in pickle_double:
        os.remove(pd[0]) 
                    
    #end = timer()
    #print(f'elapsed time (extraction des dictionnaires du pickle et combinaison): {end - start}')
    
    #for trip, FW in shared_dict.items():
    #    print("{}:   {}".format(trip, FW))
    
    return dict(shared_dict), taxanames


def del_replications(treerep, method="VR", verbose=False): #pour un arbre treerep

    def del_replications_node(treerep, method="VR", verbose=False): 
        """
            prend un arbre en argument 
            vérifie si réplication
            si oui supression de la réplication pour un seul noeud
            (dans l'ordre, apical, puis orthologue, puis paralogue)
            puis renvoi de l'arbre (apical ou orthologue) ou des 
            sous-arbres (paralogue) en deux listes: la première
            contenant les arbres avec répétitions, l'autre
            contenant les arbres sans répétitions
        """
        paralogs = dict()
        cardl = list(Tree.get_leaf_names(treerep))
        cards = set(Tree.get_leaf_names(treerep))
        cardr = [x for x in cards if cardl.count(x) > 1] #liste des taxons répétés
        treelistrep = []
        treelistnorep = []
        treelistverif = []
        
        if cardr: #détection de répétitions
            nodelist = sorted([node for node in treerep.traverse(strategy="postorder") if node.is_leaf() == False], key=lambda n: len(n))
    
            for e in nodelist: #itère dans l'ordre du plus petit noeud au plus grand (en nombre de feuilles)
                leafset = set(Tree.get_leaf_names(e))
                
                #traitement du noeud le plus petit contenant des terminaux répétés
                if len(leafset) != len(Tree.get_leaf_names(e)): #si présence de répétitions
                
                    #identifier si la répétition: apical, orthologue, paralogue
                    i = 0 #i est le nombre de noeuds internes connectés au noeud (e)
                    e.add_feature("FP", "main")
                    for child_node in e.get_children():
                        if not child_node.is_leaf():
                            i += 1
                            for cnode in child_node.traverse(strategy="postorder"):
                                cnode.add_feature("FP", i) #note tous les noeuds à sauver à l'intérieur de main
    
                    #apical
                    if i == 0:
                        paralogs["main"] = treerep.copy(method='cpickle')
                        for l in set(cardr) & leafset: #pour chaque terminal
                            leaf1 = False
                            for delnode in paralogs["main"].search_nodes(FP="main")[0].get_leaves_by_name(name=l):
                                if leaf1 == False:
                                    leaf1 = True #sauf la première                            
                                else:
                                    delnode.delete() #supprime les feuilles
                                    if verbose:
                                        print("apical repetition:{}".format(delnode.name))
    
                    #orthologue
                    elif i == 1:
                        paralogs["main"] = treerep.copy(method='cpickle')
                        kl = []
                        
                        for l in set(cardr) & leafset: #pour chaque terminal
                            
                            leaf1 = False
                            
                            for sorted_nodes in sorted([node for node in paralogs["main"].search_nodes(FP="main")[0].traverse(strategy="postorder") if node.is_leaf() == False], key=lambda n: len(n)):
                                for childnode in sorted_nodes.get_children():
                                    if childnode.is_leaf():
                                        if childnode.name == l: #pour chaque feuille répétée sauf la première
                                            if leaf1 == False:
                                                leaf1 = True
                                            elif leaf1 == True:
                                                kl.append(childnode)
            
                        for delnode in kl:
                            delnode.delete() #supprime les feuilles
                            if verbose:
                                print("ortholog repetition:{}".format(delnode.name))

                    #paralogue
                    else:
                        if method == "VR": #uniquement dans la procédure de Rineau
                            paralogs["main"] = treerep.copy(method='cpickle')
                            
                            for delnode in paralogs["main"].search_nodes(FP="main")[0].iter_descendants("postorder"):
                                if not delnode.is_leaf():
                                    delnode.delete() #supprime les noeuds internes inclus dans le paralogue ou il y a répétition
            
                        j = 1
                        while j != i+1:
                            if method == "VR":
                                paralogs[j] = e.copy(method='cpickle') #copie du noeud paralogue uniquement
                            elif method == "RZB":
                                paralogs[j] = treerep.copy(method='cpickle') #copie du noeud paralogue uniquement
    
                            for delnode in paralogs[j].search_nodes(FP="main")[0].iter_descendants("postorder"): 
                                if not delnode.is_leaf() and not delnode.FP == j:
                                    delnode.delete() #supprime les noeuds internes inclus dans le paralogue ou il y a répétition
        
                            j += 1
                        
                        if verbose:
                            print("paralog repetition:{}.{} subtrees generated from the main tree.".format(set(cardr) & leafset,i+1))

                    break #la fonction ne traite qu'une seule instance
            
            treelistverif = [p for p in paralogs.values()]
        
        else: #si pas de répétition
            treelistnorep.append(treerep)
        
        #tri entre arbres répétés et non-répétés
        for t in treelistverif:
            if len(set(Tree.get_leaf_names(t))) > 2:
                if len(list(Tree.get_leaf_names(t))) != len(set(Tree.get_leaf_names(t))): #détection de répétitions
                    treelistrep.append(t)
                else:
                    treelistnorep.append(t)
                    
        #supression des attributs (main) : copie propre des arbres pour traitement récursif
        treelistrepclearcopy = []
        
        for t in treelistrep:
            treelistrepclearcopy.append(t.copy(method="newick"))

        return treelistrepclearcopy, treelistnorep

    ########################

    listundone, listdone = del_replications_node(treerep, method, verbose)
    
    if listundone: #si répétition
        while listundone: #boucle tant que listundone contient des arbres
            listundone1 = []
            listdone1 = []    
        
            for treerep2 in listundone:
                LU, LD = del_replications_node(treerep2, method, verbose)
                listundone1 += LU
                listdone1 += LD
        
            #à la fin de cette boucle, tous les arbres de listundone ont été traités
            listundone = listundone1.copy() #on remplace les arbres avec répétition par la nouvelle liste
            listdone += listdone1 #par contre pour les arbres terminés, on les accumule dans listdone
            
    #détection des arbres non-informatifs
    treelist = []
    for infotree in listdone:
        
        #détection et supression des noeuds internes vides (un seul descendant)
        for node in infotree.traverse():
            if not node.is_leaf() and len(node.get_children()) == 1:
                node.get_children()[0].delete()
        
        #sélection de l'arbre seulement s'il a plus d'un noeud interne
        if len([t for t in infotree.traverse() if not t.is_leaf()]) > 1:
            infotree.ladderize()
            treelist.append(infotree)
    
    if verbose: #bloc à afficher si flag verbose     
    
        print("output subtrees:")
        if treelist:
            for l in treelist:
                print(l.write(format=9))
        else:
            print("No informative tree")

    return treelist #retourne une liste d'arbres sans répétition
     