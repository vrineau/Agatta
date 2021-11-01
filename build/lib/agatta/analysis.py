# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:42 2020

@author: Valentin Rineau
"""

from .ini import character_extraction, taxa_to_numbers
from fractions import Fraction
from itertools import combinations, product
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from functools import partial
from tqdm import tqdm
import warnings
import pickle
import uuid
import shutil
import os
import sys
import time

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
    from ete3 import Tree
    


class triplet():
    __slots__ = ("in_taxa", "out_taxa", "FW")

    def __init__(self, in_taxa = set(), out_taxa = set(), FW = Fraction()):
        self.in_taxa = set(in_taxa)
        self.out_taxa = set(out_taxa) # taxon out
        self.FW = FW # Liste contenant la fraction de FW, avec d'abord le numérateur puis le dénominateur
   
    def __eq__(self, other_triplet):
        return (set(self.in_taxa) == set(other_triplet.in_taxa)) and (self.out_taxa == other_triplet.out_taxa)

    def __call__(self, in_taxa = set(), out_taxa = set(), FW = Fraction(), parent_triplets = ()):
        return self

    def __hash__(self):
        return hash(list(self.in_taxa)[0]) ^ hash(list(self.in_taxa)[1]) ^ hash(list(self.out_taxa)[0]) 

    def __del__(self):
        del self

    def __repr__(self):
        return "({}({},{}))".format(str(list(self.out_taxa)[0]), 
                                    str(list(self.in_taxa)[0]), 
                                    str(list(self.in_taxa)[1]))


def triplet_extraction(infile, taxa_replacement_file=False): 
    """
        taxa replacement n'est pas encore implémenté'
    """
    
    print("Loading triplet set")
    
    triplet_dict = {} #dictionnaire ou sont enregistrés l'ensemble des arbres de caractères du fichier "pathmessage"
    
    if taxa_replacement_file:
        taxa_dict = {} #dictionnaire de conversion des taxons/symboles
        with open(taxa_replacement_file, "r") as line:
            taxa_dict[int(line.split("    ")[0])] = line.split("    ")[0]

    with open(infile, "r") as file_tree :
        for line in file_tree: 
            tripletstr = line.split("    ")
            
            newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in tripletstr[0])
            taxaint = [int(i) for i in newstr.split()]
            trip = triplet({taxaint[1], taxaint[2]}, {taxaint[0]})
            
            if taxa_replacement_file:
                convert_trip = triplet({taxa_dict[taxaint[1]], 
                                        taxa_dict[taxaint[2]]}, 
                                       {taxa_dict[taxaint[0]]})
                triplet_dict[convert_trip] = Fraction(tripletstr[1])
            else:
                triplet_dict[trip] = Fraction(tripletstr[1])
                    
    print("{} characters loaded".format(str(len(triplet_dict))))
    
    return triplet_dict


def picklemerge(namelist):
    """
    prend une liste de noms, ouvre les pickle correspondant et combine les dicos 
    contenus dedans en un superdict
    """
    with open(namelist[0], 'rb') as pickle_file:
        superdict = defaultdict(int,pickle.load(pickle_file))
                
    for i in range(1,len(namelist)):
        with open(namelist[i], 'rb') as pickle_file:
            tripletdict_temp = pickle.load(pickle_file)
            for trip, FW in tripletdict_temp.items():
                superdict[trip] += FW
                
    return superdict   
    

def tripdecFW(triplet_output=dict(), total_taxaset=False, character=Tree()):
    """
        prend un arbre et renvoie une liste de triplets / FW

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
                if not triplet({taxa_in1, taxa_in2}, {singleout}) in tree_tripdic:
                    nodaltripletdict[triplet({taxa_in1, taxa_in2}, {singleout})] = 0
                
        for trip in nodaltripletdict.keys():
            tree_tripdic[trip] += Fraction(totalnodeweight, len(nodaltripletdict))

        return tree_tripdic


def tripdec(weighting, character):
    """
        prend un arbre et renvoie une liste de triplets / FWNL MW AW NW
    """
    #calcul des triplets par caractère
    

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
                    tree_tripdic[triplet({taxa_in1, taxa_in2}, {singleout})] += FW_node
                elif weighting == "UW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2}, {singleout})] += 1
                if weighting == "MW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2}, {singleout})] = MW_node
                elif weighting == "AW" or weighting == "NW":
                    tree_tripdic[triplet({taxa_in1, taxa_in2}, {singleout})] = 1

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


def tripdec_allweights(weighting, character): #CHANGER DE NOM
    """
        prend un arbre et renvoie une liste de triplets / FW FWNL MW AW NW
        Appelle tripdec et tripdecFW
        Utilisé uniquement par parallel tripdec
    """
    #calcul des triplets par caractère

    pickle_dict = defaultdict(int)
    tree_id = str(uuid.uuid4())
        
    cardinal_character1 = set(Tree.get_leaf_names(character)) #obtenir la liste des taxons de l'arbre
    internal_nodes1 = (node for node in character.traverse(strategy="postorder") if node.is_leaf() == False and node.is_root() == False)
    card = len(cardinal_character1)
    MW_node = Fraction(6,card*(card-1))
    #calculer les triplets par noeud => construire une liste de triplets
        
    for node in internal_nodes1:
        pickle_name = str(uuid.uuid4()) + "_agatta.pickle"
        nodaltripletdict = dict()
        taxa_in = set(Tree.get_leaf_names(node)) #obtenir la liste des taxons d'un noeud - IN
        taxa_out = cardinal_character1 - taxa_in #obtenir la liste des taxons hors d'un noeud - OUT  
        # print(str(len(cardinal_character1)))
        # print(str(len(taxa_in)))
        # print(str(len(taxa_out)))
        # print(str(node.is_root()))
        FWNL_node = Fraction(2, len(taxa_in))
        tempw = sum((len(t1)*len(t2) for t1, t2 in combinations((child_node.get_leaf_names() for child_node in node.get_children()), 2)))
        FW_node = Fraction(len(taxa_out)*(len(taxa_in)-1), len(taxa_out)*tempw)

        for singleout in taxa_out: #construction des triplets et ajout au set du noeud  
            
            if weighting in ("FW","MW","AW","NW"):
                
                for taxalist1, taxalist2 in combinations((child_node.get_leaf_names() for child_node in node.get_children()), 2):
                    for taxa_in1, taxa_in2 in product(taxalist1, taxalist2):
                                                
                        if weighting == "FW":
                            nodaltripletdict[triplet(in_taxa={taxa_in1,taxa_in2},out_taxa={singleout})] = FW_node
            
                        elif weighting == "MW":
                            nodaltripletdict[triplet(in_taxa={taxa_in1,taxa_in2},out_taxa={singleout})] = MW_node
                            
                        else: #AW NW
                            nodaltripletdict[triplet(in_taxa={taxa_in1,taxa_in2},out_taxa={singleout})] = 1            

            else: #FWNL ou UW
                    
                for taxa_in1, taxa_in2 in combinations(taxa_in, 2): #obtenir les combinaisons in_taxa
                                 
                    if weighting == "FWNL":
                        nodaltripletdict[triplet(in_taxa={taxa_in1,taxa_in2},out_taxa={singleout})] = FWNL_node
                        
                    else: # UW
                        nodaltripletdict[triplet(in_taxa={taxa_in1,taxa_in2},out_taxa={singleout})] = 1
                    

        with open(pickle_name, 'wb') as pickle_file:
            pickle.dump(nodaltripletdict, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
            
        pickle_dict[pickle_name] = tree_id

    return pickle_dict


def standard_tripdec(character_dict, weighting, prefix=False, verbose=True):
    """
    Prend un tuple d'arbres et compute toutes les listes de triplets en parallèle
    Combine les listes en une seule
    Renvoie la liste
    
    Fonctionne pour tous types de pondérations

    """
    

    #replace taxa by numbers
    taxa_dict, code_dict, taxa_convers = taxa_to_numbers(character_dict)
    
    triplet_output = dict() 
    
    if weighting == "FW":
        
        for treedec in character_dict.keys():
            triplet_output2 = tripdecFW(triplet_output, False, character=treedec)
            
            for trip, FW in triplet_output2.items():
                if trip in triplet_output: #nouveaux triplets
                    triplet_output[trip] += FW  
                else: #incrémentation de triplets préexistants
                    triplet_output[trip] = FW
                
    else: #MW, AW, NW
        for treedec in character_dict.keys():
            triplet_output2 = tripdec(weighting, treedec)
            
            for trip, FW in triplet_output2.items():
                if trip in triplet_output: #nouveaux triplets
                    triplet_output[trip] += FW  
                else: #incrémentation de triplets préexistants
                    triplet_output[trip] = FW

    if prefix:
        with open(prefix+".triplet", "w") as tdfile:
            for trip, FW in triplet_output.items():
                
                #remettre les noms str
                trip2 = trip
                in_taxa = list(trip2.in_taxa)
                trip2.in_taxa = {taxa_dict[in_taxa[0]], taxa_dict[in_taxa[1]]}
                trip2.out_taxa = {taxa_dict[list(trip.out_taxa)[0]]}
                
                tdfile.write("{}:    {}    {}\n".format(trip, FW, round(float(FW),4)))        


        with open(prefix+".taxabloc", "w") as taxa_bloc_file:
            for taxa, code in taxa_dict.items():
                taxa_bloc_file.write("{}    {}\n".format(str(code),taxa))
                
    return triplet_output


def parallel_tripdec(character_dict, weighting, prefix=False, ncpu="auto"):
    """
    Prend un tuple d'arbres et compute toutes les listes de triplets en parallèle
    Combine les listes en une seule
    Renvoie la liste
    
    Fonctionne pour tous types de pondérations

    """
    
    def split(a, n):
        """
        prend une liste de pickle names et les trie en n listes
        """
        k, m = divmod(len(a), n)
        return list(a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


    #nombre de coeurs    
    if ncpu == "auto":
        ncpu = cpu_count()
    else:
        ncpu = int(ncpu)
        
    os.mkdir("out_agatta")
    os.chdir("out_agatta")

    #computation des dictionnaires de triplets pour chaque arbre en parallèle
    treetuple = tuple(character_dict) #dictionnaire en tuple pour multiprocessing
    
    # #création du dictionnaire 
    # taxanames = list(set([x for l in [character.get_leaf_names() for character in character_dict] for x in l]))
    
    with Pool(ncpu) as pool: #activation du multiprocessing
        func = partial(tripdec_allweights, weighting)
        pickle_double = pool.map(func, treetuple) #liste de dictionnaires de triplets

    #étape d'extraction des dictionnaires du pickle et combinaison    
    picknum = len([item for sublist in [tuple(dicname.keys()) for dicname in pickle_double] for item in sublist])
    
    if ncpu > picknum: #if more cpu than pickles, use the minimal number of cpu
        ncpu = picknum
        pickletuple = split([item for sublist in [tuple(dicname.keys()) for dicname in pickle_double] for item in sublist], ncpu)
    else: #else, use the maximum number of cpu allowed
        pickletuple = split([item for sublist in [tuple(dicname.keys()) for dicname in pickle_double] for item in sublist], ncpu)
        
    with Pool(ncpu) as pool: #activation du multiprocessing
        dictuple = pool.map(picklemerge, pickletuple) #liste de dictionnaires de triplets
 
    shared_dict = dictuple[0]
        
    if weighting == "NW": #poids unique de 1

        for d in dictuple[1:len(dictuple)]:
            shared_dict.update(d)
 
    else: #somme des poids

        for d in dictuple[1:len(dictuple)]:
            for trip, FW in d.items():
                shared_dict[trip] += FW
                         
    #supression des pickle
    os.chdir('..')
    shutil.rmtree(os.getcwd()+'/out_agatta', ignore_errors=True)

    if prefix:
        with open(prefix+".triplet", "w") as tdfile:
            
            #replace taxa by numbers
            taxa_dict, code_dict, taxa_convers = taxa_to_numbers(character_dict)
            
            for trip, FW in shared_dict.items():
            
                in_taxa = list(trip.in_taxa)

                tdfile.write("({},({},{})):    {}    {}\n".format(list(trip.out_taxa)[0], 
                                                                  in_taxa[0], in_taxa[0],
                                                                  FW, round(float(FW),4)))        


        with open(prefix+".taxabloc", "w") as taxa_bloc_file:
            for taxa, code in taxa_dict.items():
                taxa_bloc_file.write("{}    {}\n".format(str(code),taxa))

    return shared_dict


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
        
        sys.setrecursionlimit(10000) #ajout d'un except?
        
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
                                
                    for node in treerep.traverse():
                        try:
                            if node.FP:
                                pass
                        except AttributeError:
                            node.add_feature("FP", "out")

                    #cas spécial: repeated leaf branched to a symmetric node
                    paralog_but_leaf = False
                    
                    if i > 1:                    
                        special_tree = treerep.copy(method='cpickle')
                        
                        for child in special_tree.search_nodes(FP="main")[0].get_children():
                            if child.is_leaf:
                                if child.name in cardr:
                                    paralog_but_leaf = True
                                    print("ortholog repetition (on symmetric node):{}".format(child.name))
                                    child.delete()
                        
                    #cas spécial feuilles sur paralogue
                    if paralog_but_leaf:
                        paralogs["main"] = special_tree
                    
                    #apical
                    elif i == 0:
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

                        j = 1
                        
                        #procedure Rineau
                        if method == "VR": #uniquement dans la procédure de Rineau
                            paralogs["main"] = treerep.copy(method='cpickle')
                            
                            for delnode in paralogs["main"].search_nodes(FP="main")[0].iter_descendants("postorder"):
                                if not delnode.is_leaf():
                                    delnode.delete() #supprime les noeuds internes inclus dans le paralogue ou il y a répétition
            
                            while j != i+1:
                                paralogs[j] = treerep.copy(method='cpickle') #e.copy(method='cpickle') #copie du noeud paralogue uniquement
                                
                                for delnode in paralogs[j].traverse(): 
                                    if not delnode.is_leaf() and not delnode.FP == j:
                                        delnode.delete() #supprime les noeuds internes inclus dans le paralogue ou il y a répétition
            
                                j += 1

                        #procedure Nelson / Zaragueta
                        elif method == "RZB": #uniquement dans la procédure de Rineau
                            
                            while j != i+1: #pour chaque noeud branché au paralogue
                            
                                #construction d'un arbre pour chaque noeud branché au paralogue
                                paralogs[j] = treerep.copy(method='cpickle') #copie de tout
                                    
                                for delnode in paralogs[j].search_nodes(FP="main")[0].get_children(): 
                                    if not delnode.FP == j: #on sauve keepnode et detach le reste
                                        delnode.detach() #supprime les noeuds internes inclus dans le paralogue ou il y a répétition

                                paralogs[j].search_nodes(FP="main")[0].delete()
                                j += 1
                        
                        if verbose:
                            print("paralog repetition. {} subtrees generated from the main tree.".format(str(len(paralogs))))

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
            
    return treelist

def del_replications_forest(character_dict, method="VR", prefix="agatta_del_replications", verbose=False):
    
    print("Managing polymorphism")
    
    tree_dict = dict()
    
    output_trees = 0
    
    with open(prefix+".poly", "w") as logfile:
        for treed, treeid in tqdm(character_dict.items()):
            treelist = del_replications(treed, method, verbose)
            
            if treelist:
                if len(treelist) == 1:
                    tree_dict[treelist[0]] = treeid
                    logfile.write("["+str(treeid)+"] "+treelist[0].write(format=9)+"\n")
                else:
                    i = 1
                    output_trees += len(treelist)
                    for treel in treelist:
                        tree_dict[treel] = str(treeid)+"."+str(i)
                        
                        logfile.write("["+str(treeid)+"."+str(i)+"] "+treel.write(format=9)+"\n")
                        i+= 1
            else:
                logfile.write("["+str(treeid)+"] no informative tree\n")         

    if output_trees != 0:
        print("Polymorphism removed. "
              "{} polymorphism free informative trees computed.".format(
                  str(output_trees)))
    else:
        print("Polymorphism removed. No polymorphism free "
              "informative tree remaining.")

    

    return tree_dict #retourne une liste d'arbres sans répétition


def rep_detector(character_dict):
    
    for t in character_dict.keys():
        cardl = list(Tree.get_leaf_names(t))
        cards = set(Tree.get_leaf_names(t))
        
        if len(cardl) != len(cards):
            return True
    
    return False
    

def main_tripdec(input_item, prefix, taxa_replacement, weighting, parallel, verbose):
    """
    

    Parameters
    ----------
    input_item : string or dictionary
        can be the path of the file or the dictionary of trees.
    prefix : TYPE
        DESCRIPTION.
    taxa_replacement : TYPE
        DESCRIPTION.
    weighting : TYPE
        DESCRIPTION.
    parallel : TYPE
        DESCRIPTION.
    verbose : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    if type(input_item) == str:
        input_item = character_extraction(input_item, taxa_replacement=taxa_replacement)
    
    #abort if repetitions
    if rep_detector(input_item):
        raise UserWarning('Characters contain repetitions. Operation aborted. Please remove repetitions before conversion.')

    #computation des dictionnaires de triplets pour chaque arbre en parallèle
    if verbose:
        print("Starting triplet decomposition and {} computation".format(weighting))
        start = time.time()
    
    if parallel == "no":
        triplet_dict = standard_tripdec(input_item, weighting, prefix)
    else:
        triplet_dict = parallel_tripdec(input_item, weighting, prefix, parallel)
        
    if verbose:    
        end = time.time()
        time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
        print("elapsed time (triplet decomposition and weighting): {}".format(time_cptr))
        print(str(len(triplet_dict))+" triplets decomposed")
    

    return triplet_dict