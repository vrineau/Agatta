# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:52 2020

@author: Valentin Rineau
"""

from fractions import Fraction
from itertools import combinations
from tkinter import filedialog
from tkinter import Tk
from .analysis import triplet, standard_tripdec
from .ini import taxa_list_extraction
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
    from ete3 import Tree


def constrict(treelist, prefix=False): 
    """
        prend une liste d'arbres et retourne un dictionnaire contenant 
        le consensus strict uniquement
    """
    
    print("Strict consensus computation")
    
    constrict = treelist[0] #arbre qui servira de base pour le consensus strict (peu importe lequel)
    del treelist[0] #on supprime le premier arbre qui sert de toute façon de base
    kill_node_list = [] #liste des noeuds à supprimer à la toute fin de l'analyse de consensus
    
    for node in constrict.traverse(strategy="preorder"): #pour chaque noeud informatif de l'arbre...
        if not node.is_leaf() and not node.is_root():  #...on va devoir vérifier qu'il existe chez chacun des autres arbres
            for tree_compare in treelist: #pour chaque arbre de la liste...
                node_find = False #compteur du noeud à vérifier
                for node_tree_compare in tree_compare.traverse(strategy="preorder"): #...on regarde s'il y a le noeud
                    if set([leaf.name for leaf in node.get_leaves()]) == set([leaf.name for leaf in node_tree_compare.get_leaves()]):
                        node_find = True
                        break #on peut passer à l'arbre suivant
                if node_find == False: #à la fin de l'analyse de chaque arbre, si le noeud n'a pas été retrouvé, on le supprime chez le 
                    kill_node_list.append(node)
                    break

    for kill_node in kill_node_list: #on supprime les neuds hors de la boucle précédente pour éviter de modifier l'itération elle-même
        kill_node.delete()
    
    if prefix:
        with open(prefix+".constrict", "w") as constrictfile:
            constrictfile.write(constrict.write(format=9)+"\n")

    print("Strict consensus computed")
    
    return constrict

def rcc(treelist, prefix=False, verbose=False): 
    """
        Prend un dictionnaire d'arbres et renvoie une liste contenant 
        tous les arbres du profil du reduced cladistic consensus 
        (RCC; Wilkinson 1994)
    """
    
    def nts_to_trees(cardinal, in_list): 
        """
            Arguments: cardinal & sets in_taxa en value
        """
    
        t = Tree() #créer un arbre
        t.add_features(taxa_content = cardinal)
        
        for nts in sorted(in_list, key=lambda x: len(x), reverse=True): #on ajoute du plus grand au plus petit nts
            for node_iter in t.traverse(strategy="postorder"): #rechercher tous les noeuds, du moins inclusif au plus inclusif
                
                #récupérer le set des taxons de tous les noeuds fils directs
                taxa_set = set() 
                for node_child1 in node_iter.get_children(): 
                    taxa_set = taxa_set | node_child1.taxa_content
                                    
                #si il y a de la place pour ajouter un noeud fils => vérification de tous les nts par taille => si match, ajout de noeud fils
                if len(node_iter.taxa_content - taxa_set) > 1 and ((node_iter.taxa_content - taxa_set) >= nts or (node_iter.taxa_content - taxa_set) == nts): #si encore possible qu'il manque un noeud fils à node_iter
                    node_iter.add_child().add_features(taxa_content = nts) #créer un noeud fils avec in_taxa
                    break
        
        #label des noeuds: orthologue, paralogue, apical
        for node in t.traverse(strategy="preorder"):
            if [node1 for node1 in node.get_children() if node1.taxa_content]: #s'il y a des noeuds fils
                #récupérer le set des taxons de tous les noeuds fils directs
                taxa_set = set() 
                for node_child2 in node.get_children(): 
                    taxa_set = taxa_set | node_child2.taxa_content
                    
                for taxa in node.taxa_content - taxa_set:
                    node.add_child(name=taxa).add_features(taxa_content = set()) #ajout d'une feuille pour chaque membre du set
                    
            elif node.taxa_content: #s'il n'y a pas de noeuds fils
                for taxa in node.taxa_content:
                    node.add_child(name=taxa).add_features(taxa_content = set()) #ajout d'une feuille pour chaque membre du set
        
        return t #fin de la fonction, output = arbre


    #fonction de tri des arbres
    def get_len(in_set):
        return len(in_set)

    print("Reduced cladistic consensus (RCC) computation")

    #vérification de la taille des arbres
    cardinal_list = []
    for treel in treelist:
        leaf_cptr = 0
        for node in treel.traverse(strategy="preorder"):
            if node.is_leaf() == True:
                leaf_cptr += 1
        cardinal_list.append(leaf_cptr)
        
    cardinal_set = set(cardinal_list)
    
    if len(cardinal_set) > 1:
        raise ValueError("les arbres n'ont pas tous le même nombre de terminaux")
    
    #construction de la liste de composantes par arbre
    component_trees = []
        
    for treel in treelist:
        component_trees.append([]) #nouvelle liste correspondant à la liste de composantes d'un arbre
        for node in treel.traverse(strategy="preorder"):
            if not node.is_leaf() and not node.is_root():
                component_trees[-1].append([set([leaf.name for leaf in node.get_leaves()]),set([leaf.name for leaf in treel.get_leaves()]) - set([leaf.name for leaf in node.get_leaves()])]) #ajout des composantes de structure (in,out)
    
    #construction de la liste d'intersection *nts_intersect*
    nts_intersect = []
    
    for comp1 in component_trees[0]:
        for comp2 in component_trees[1]:
            if len(comp1[0] & comp2[0]) > 1 and len(comp1[1] & comp2[1]) > 0:
                nts_intersect.append([comp1[0] & comp2[0],comp1[1] & comp2[1]])
    
    del component_trees[1] #supression des arbres de component_trees après transfert vers nts_intersect
    del component_trees[0] #supression des arbres de component_trees après transfert vers nts_intersect
    
    #supression de la redondance
    del_nts_intersect = [] #liste ou sont stockés les nts à suprimer
    for nts1, nts2 in combinations(nts_intersect, 2): #obtenir les combinaisons de nts de la liste d'intersection
        if (nts1[0] <= nts2[0] and (nts1[1] <= nts2[1] or nts1[1] == nts2[1])) and (nts1 not in del_nts_intersect): #si nts1 est un sous-arbre de nts2 et que nts1 n'est pas déjà dans la liste de supression
            del_nts_intersect.append(nts1) #ajouter nts1 à la liste des nts à supprimer
            
        elif (nts1[0] >= nts2[0] and (nts1[1] >= nts2[1] or nts1[1] == nts2[1])) and (nts2 not in del_nts_intersect): #si nts2 est un sous-arbre de nts1 et que nts1 n'est pas déjà dans la liste de supression
            del_nts_intersect.append(nts2) #ajouter nts1 à la liste des nts à supprimer
    
    for del_nts in del_nts_intersect: #boucle de supression des nts redondants détectés
        while del_nts in nts_intersect:
            nts_intersect.remove(del_nts)
    
    #intersection de chaque ensemble de composantes (par arbre) avec la liste d'intersection *nts_intersect*
    while len(component_trees) > 0:
        nts_intersect_new = []
        for comp1 in component_trees[-1]: #algorithme à pile
            for comp2 in nts_intersect:
                if len(comp1[0] & comp2[0]) > 1 and len(comp1[1] & comp2[1]) > 0:
                    nts_intersect_new.append([comp1[0] & comp2[0],comp1[1] & comp2[1]])
    
        nts_intersect = nts_intersect_new
    
        del component_trees[-1] #supression du pool de composantes utilisé
        
        #supression de la redondance
        del_nts_intersect = [] #liste ou sont stockés les nts à suprimer
        doubles_nts_intersect = []
        for nts1, nts2 in combinations(nts_intersect, 2): #obtenir les combinaisons de nts de la liste d'intersection
            if nts1 == nts2:
                doubles_nts_intersect.append(nts1)
                
            elif (nts1[0] <= nts2[0] and (nts1[1] <= nts2[1] or nts1[1] == nts2[1])) and (nts1 not in del_nts_intersect): #si nts1 est un sous-arbre de nts2 et que nts1 n'est pas déjà dans la liste de supression
                del_nts_intersect.append(nts1) #ajouter nts1 à la liste des nts à supprimer
                
            elif (nts1[0] >= nts2[0] and (nts1[1] >= nts2[1] or nts1[1] == nts2[1])) and (nts2 not in del_nts_intersect): #si nts2 est un sous-arbre de nts1 et que nts1 n'est pas déjà dans la liste de supression
                del_nts_intersect.append(nts2) #ajouter nts1 à la liste des nts à supprimer
                
        for del_nts in doubles_nts_intersect: #boucle de supression des nts redondants détectés
            while nts_intersect.count(del_nts) != 1:
                nts_intersect.remove(del_nts)
    
        for del_nts in del_nts_intersect: #boucle de supression des sous-arbres
            while del_nts in nts_intersect:
                nts_intersect.remove(del_nts)
    
        if len(nts_intersect) == 0: #si nts_intersect est vide, le résultat sera un rateau, on peut gagner du temps et break le while
            break
    
    #regrouper les nts en liste de cardinal commun
    nts_groups = {}
    for nts in nts_intersect: #ajout de l'information du cardinal en position [2] de chaque nts (liste)
        nts.append(frozenset(nts[0]|nts[1]))
        nts_groups[frozenset(nts[0]|nts[1])] = []
    
    for nts in nts_intersect: #groupement des nts de même cardinal
        nts_groups[nts[2]].append(nts[0]) #ajouter les taxons in
    
    for nts in nts_intersect: #ajout du maximum d'information possible (intersection des nts)
        for cardinal, in_list in nts_groups.items():
            if nts[2] != cardinal and nts[2] >= cardinal and len(nts[0] & set(cardinal)) > 1:
                nts_groups[cardinal].append(nts[0] & set(cardinal))
    
    #fusion des nts en arbres
    profile = []
    
    for cardinal, in_list in nts_groups.items(): #pour chaque arbre (chaque ensemble de noeuds)
        tree_profile = nts_to_trees(cardinal, in_list)
        if not len(tree_profile) == leaf_cptr:
            profile.append(tree_profile)
        
    profile.append(constrict(treelist, prefix=False, verbose=False))
    
    #tri des arbres par taille
    profile = sorted(profile, key=get_len)
    
    if verbose:
        if not cardinal_list[0] in [len(profile_tree) for profile_tree in profile]:
            print("Le profil RCC ne contient pas le consensus strict (i.e. le consensus strict n'est pas informatif)")
            
        #affichage des résultats
        if len(nts_intersect) == 0:
            print("Profil RCC non-informatif (il n'y a pas d'information phylogénétique commune entre les arbres)")
        else:
            for tree_profile in profile:
                print(tree_profile)
            
    if prefix:
        with open(prefix+".rcc", "w") as rccfile:
            first_line = True
            i = 1
            for profiletree in profile:
                if first_line:
                    rccfile.write("Strict consensus:    "+profiletree.write(format=9)+"\n")
                    first_line = False
                else:
                    rccfile.write("Profile tree "+str(i)+":    "+profiletree.write(format=9)+"\n")
                    i += 1

    print("RCC of {} trees computed".format(str(len(profile))))
    
    return profile

def RI(cladogram_dict, character_dict, weighting="FW", output=None):
    """
        Fonction permettant de calculer l'indice de rétention global (total)
        et pour chaque caractère.
        prend en argument un dictionnaire de caractères et 
        un dictionnaire cladogramme
    """
    
    print("Computing retention index")
    
    c_triplet_dict = standard_tripdec(cladogram_dict, 
                                      weighting,
                                      prefix=False,
                                      verbose=False)

    RI_tot = [0, 0] #numérateur, dénominateur de l'IR total
    RI_char_dict = {} #IR par arbre
    RI_char_dict_num = {} #numérateur par arbre
    RI_char_dict_denom = {} #dénominateur par arbre

    #calcul de l'IR par caractère et total en même temps
    for chartree, keys in character_dict.items():
        
        RI_char_dict[keys] = 0
        RI_char_dict_num[keys] = 0
        RI_char_dict_denom[keys] = 0
        
        triplet_dict = standard_tripdec({chartree: keys},
                                        weighting,
                                        prefix=False,
                                        verbose=False)
        
        for triplet1, FW in triplet_dict.items(): #pour chaque triplet à l'intérieur d'un caractère
            #IR numérateur
            if triplet1 in c_triplet_dict: #si le triplet du caractère est dans le cladogramme, on ajoute la valeur à l'IR
                
                RI_char_dict_num[keys] += FW #par caractère
                RI_tot[0] += FW #pour l'IR total
                
            #IR dénominateur
            RI_char_dict_denom[keys] += FW
            RI_tot[1] += FW # construction du dénominateur de l'IR total
            
    for keys, values in RI_char_dict_num.items():
        
        if RI_char_dict_denom[keys] == 0:
            RI_char_dict[keys] = Fraction(0, 1) 
            
        else:
            RI_char_dict[keys] = Fraction(RI_char_dict_num[keys], RI_char_dict_denom[keys]) #construction des IR par ajout du dénominateur
    
    if RI_tot[1] == 0:
        RI_char_dict["Total retention index"] = Fraction(0, 1) # ajout de l'IR  total au dictionnaire
    else:
        RI_char_dict["Total retention index"] = Fraction(RI_tot[0], RI_tot[1]) # ajout de l'IR  total au dictionnaire

    print("Retention index computed")
    
    #génération du fichier nexus
    if output == None:
        root = Tk()
        root.filename =  filedialog.asksaveasfilename(title = "Save retention index file", filetypes = (("text files","*.txt"),("all files","*.*")))
        root.withdraw()

    with open(output+".ri", "w") as RI_file:
        RI_file.write("#Retention index\n")
    
        for keys, values in RI_char_dict.items(): #affichage de la liste des IR
            RI_string = str(keys)+": "+str(round((float(values)*100), 2))
            print(RI_string)
            RI_file.write("\n"+RI_string)
            
            
def ITRI(true_tree, reconstructed_tree, prefix, weighting="FW"):  
    """
        prend en argument deux dictionnaires de UN caractère (un pour chacun
        des deux arbres dont on souhaite calculer la distance ITRI)
        liste de triplets et un dictionnaire cladogramme/liste de triplets => fonction
        triplet_decomposition_d(character_extraction())
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
    RI_intersect_reconstructed_tree = Fraction(0, 1) #total reconstruted tree intersection total true tree
    RI_intersect_true_tree = Fraction(0, 1) #total reconstruted tree intersection total true tree

    for trip, FW in rt_tripdic.items(): #pour chaque ensemble de triplets correspondant à un caractère
        if trip in tt_tripdic: #si le triplet du caractère est dans le cladogramme, on ajoute la valeur à l'IR
            RI_intersect_reconstructed_tree += FW #par caractère
        RI_reconstructed_tree += FW # construction du dénominateur de l'IR total
            
    for trip, FW in tt_tripdic.items(): #pour chaque ensemble de triplets correspondant à un caractère
        if trip in rt_tripdic: #si le triplet du caractère est dans le cladogramme, on ajoute la valeur à l'IR
            RI_intersect_true_tree += FW #par caractère
        RI_true_tree += FW # construction du dénominateur de l'IR total

    ITRI_power = Fraction(RI_intersect_true_tree, RI_true_tree)
    
    try:
        ITRI_arte = 1 - Fraction(RI_intersect_reconstructed_tree, RI_reconstructed_tree) #il peut y avoir l'exception ZeroDivisionError qui se lève
    except ZeroDivisionError:
        ITRI_arte  = 0
    
    ITRI_efficiency = ((float(ITRI_power)*100) - (float(ITRI_arte)*100))

    with open(prefix+".txt", "w") as itrifile:
        itrifile.write("power : "+str(round(float(ITRI_power),3))+"\n")
        itrifile.write("artefact : "+str(round(float(ITRI_arte),3))+"\n")
        itrifile.write("efficiency : "+str(round(float(ITRI_efficiency),3))+"\n")

    print("power : "+str(round(float(ITRI_power),3))+"\n")
    print("artefact : "+str(round(float(ITRI_arte),3))+"\n")
    print("efficiency : "+str(round(float(ITRI_efficiency),3))+"\n")
    
    return float(ITRI_power)*100, float(ITRI_arte)*100, ITRI_efficiency

def triplet_distance(true_tree, reconstructed_tree, prefix, method="itrisym_sum", weighting="FW"):
    """
        prend en argument deux arbres dans deux dictionnaires et calcule la triplet distance
        Contrairement à la finction ITRI, celle-ci créé une mesure de distance (symétrique)
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
    RI_intersect_reconstructed_tree = Fraction(0, 1) #total reconstruted tree intersection total true tree
    RI_intersect_true_tree = Fraction(0, 1) #total reconstruted tree intersection total true tree

    for trip, FW in rt_tripdic.items(): #pour chaque ensemble de triplets correspondant à un caractère
        if trip in tt_tripdic: #si le triplet du caractère est dans le cladogramme, on ajoute la valeur à l'IR
            RI_intersect_reconstructed_tree += FW #par caractère
        RI_reconstructed_tree += FW # construction du dénominateur de l'IR total
            
    for trip, FW in tt_tripdic.items(): #pour chaque ensemble de triplets correspondant à un caractère
        if trip in rt_tripdic: #si le triplet du caractère est dans le cladogramme, on ajoute la valeur à l'IR
            RI_intersect_true_tree += FW #par caractère
        RI_true_tree += FW # construction du dénominateur de l'IR total

    ITRIp_12 = Fraction(RI_intersect_true_tree, RI_true_tree)
    ITRIp_21 = Fraction(RI_true_tree, RI_intersect_true_tree)
    
    try:
        ITRIa_12 = 1 - Fraction(RI_intersect_reconstructed_tree, RI_reconstructed_tree) #il peut y avoir l'exception ZeroDivisionError qui se lève
    except ZeroDivisionError:
        ITRIa_12  = 0
    
    try:
        ITRIa_21 = 1 - Fraction(RI_reconstructed_tree, RI_intersect_reconstructed_tree) #il peut y avoir l'exception ZeroDivisionError qui se lève
    except ZeroDivisionError:
        ITRIa_21  = 0    
    
    ITRIe_12 = ((float(ITRIp_12)*100) - (float(ITRIa_12)*100))
    ITRIe_21 = ((float(ITRIp_21)*100) - (float(ITRIa_21)*100))

    if method == "itrisym_sum": #méthode standard pour calculer un ITRI symétrique
        ITRIsym = (ITRIe_12+ITRIe_21)/2
    
    elif method == "itrisym_product": #proposition de Grand et al. 2014
        ITRIsym = ITRIe_12*ITRIe_21
        
    with open(prefix+".txt", "w") as itrifile:
        itrifile.write(str(round(ITRIsym,3)))
    
    print("Symmetrical ITRI value :"+str(round(ITRIsym,3)))
    
    return ITRIsym
    

def character_states_test(cladogram_dict, character_dict, prefix, pdf_files=False):  
    """
        prend en argument un dictionnaire de cladogrammes et un dictionnaire de caractères
        calcule pour chaque état caractère le test proposé par Cao 2007 et
        modifié par Rineau 2017
        Si on génère les pdf, alors renseigner chemin
    """
    #import couille de ete3
    if pdf_files:
        try:
            from ete3 import NodeStyle, TreeStyle, faces, TextFace
            
            #styles de noeuds (standard, synapomorphie, homoplasie)           
            nstyle_other = NodeStyle() #style neutre
            nstyle_other["shape"] = "circle"
            nstyle_other["size"] = 0.3
            nstyle_other["fgcolor"] = "black"
            
            nstyle_accept = NodeStyle() #style synapomorphie
            nstyle_accept["shape"] = "sphere"
            nstyle_accept["size"] = 3
            nstyle_accept["fgcolor"] = "limegreen"         
            nstyle_accept.opacity = 0.3   
            
            nstyle_reject = NodeStyle() #style homoplasie
            nstyle_reject["shape"] = "sphere"
            nstyle_reject["size"] = 3
            nstyle_reject["fgcolor"] = "crimson"
            nstyle_reject.opacity = 0.3
            
            #styles de feuilles ajoutés dynamiquement (le layout est attaché au style de l'arbre)
            def lstyle(node_style):       
                #assignation des styles de feuilles                
                if node_style.is_leaf():
                    if node_style.taxa_in_state == "in":
                        name_face = TextFace(node_style.name, fgcolor="steelblue", fsize=2, bold=True)
                        #name_face.inner_background.color = "burlywood"
                    elif node_style.taxa_in_state == "out":
                        name_face = TextFace(node_style.name, fgcolor="dimgray", fsize=2)
                    else:
                        name_face = TextFace(node_style.name+" (?)", fgcolor="lightgray", fsize=2)
                    name_face.margin_left = 2
                    faces.add_face_to_node(name_face, node_style, column=0, position='branch-right')

        except ImportError:
            print("Installing PyQt5 is requested to use this functionality")
 
    print("Initiating character state test procedure")
    
    #label des noeuds du cladogramme (type et numéro, réglable avec options)
    def cladogram_label(cladogram, clade_number_option, clade_type_option):
        cladogram.ladderize()
        
        if clade_number_option == "yes": #label des noeuds du cladogramme: indispensable pour le test des états de caractères et pour la génération du pdf et du txt
            node_style_num_count = 1 
            for node in cladogram.traverse(strategy="preorder"): 
                if node.is_leaf() == False  and node.is_root() == False:
                    if pdf_files:
                        style_num = TextFace(node_style_num_count, fsize=2)
                        node.add_face(style_num, column=1, position = "branch-bottom")
                    node.add_feature("clade_label", node_style_num_count)
                    node_style_num_count += 1
                elif node.is_root() == True:
                    node.add_feature("clade_label", "root")
                else:
                    node.add_feature("clade_label", node.name)
                    
        if clade_type_option == "yes": #label des types de noeuds du cladogramme: indispensable pour le test des états de caractères
            for node in cladogram.traverse(strategy="postorder"): #label des noeuds du cladogramme
                if node.is_leaf(): #si le noeud est une feuille (ete3 fonction), on ajoute l'attribut dans clade type
                    node.add_feature("clade_type", "leaf") 
                else: #si le noeud est un noeud interne
                    i = len([child_node for child_node in node.get_children() if not child_node.is_leaf()])
                    if i >= 2:
                        node.add_feature("clade_type", "paralog")
                    elif i == 1:
                        node.add_feature("clade_type", "ortholog")
                    elif i == 0:
                        node.add_feature("clade_type", "apical_node")
                        
        syn_dict = {}
        for node_number in cladogram.traverse(strategy="postorder"):
            syn_dict[str(node_number.clade_label)] = {"accepted" : [],"rejected" : []}
                        
        return cladogram, syn_dict

    #pour une feuille, remonte les ancêtres et trouve la synapomorphie correspondante
    def find_synapomorphy(leaf_syn):
        if leaf_syn.is_root():
            ancestor_syn = leaf_syn
        elif leaf_syn.get_ancestors()[0].state == "y":
            ancestor_syn = find_synapomorphy(leaf_syn.get_ancestors()[0]) #la fonction récursive permet de boucler tant qu'on n'a pas trouvé l'ancêtre le plus inclusif
        else: #quand l'ancêtre est trouvé, il est capturé ici
            ancestor_syn = leaf_syn
        return ancestor_syn

    def synapomorphy_test_placement(cladogram, taxa_in_state, results_test_dict, char_names, char_states_names):
        results_test_dict[char_names][char_states_names] = []
        synapomorphies_set = set() #ensemble de noeuds caractérisés par l'état dérivé
        synapomorphy_test = str("accepted") #accepté/rejeté (accepté par défaut)
        potential_synapo = cladogram.get_common_ancestor(taxa_in_state) #plus petit noeud contenant tous les taxons ayant l'état dérivé
    
        #assignation de l'état au noeud:
        # si l'état est une synapomorphie sans incongruence (RI=1) ou s'il n'y a pas de paralogues impliqués
        if potential_synapo.get_leaf_names() == taxa_in_state or not potential_synapo.search_nodes(clade_type="paralog"): 
            synapomorphies_set.add(potential_synapo)
            
        #si interprétation ambigue, il faut vérifier l'arrangement des taxons au sein des paralogues
        else:
            for node_cladogram in cladogram.traverse(strategy="postorder"): #par défaut, tout le monde à state = n (y: état dérivé ou inclus dans état dérivé, sinon n)
                node_cladogram.add_feature("state", "n")
                
            for leaf_node1 in cladogram.iter_leaves(): #ajout des y à toutes les feuilles porteurs de l'état dérivé
                if leaf_node1.name in taxa_in_state:
                    leaf_node1.state = "y"
                
            for node_cladogram in cladogram.traverse(strategy="postorder"): #assignation des y aux noeuds internes (node_cladogram)
    
                if node_cladogram.clade_type == "apical_node": #si noeud apical
                    if len([apical_node_leaf for apical_node_leaf in node_cladogram.children if apical_node_leaf.state == "y"]) > 0:
                        node_cladogram.state = "y" #ajout de y si le noeud apical contient au moins une feuille y
    
                elif node_cladogram.clade_type == "paralog": #si paralogue
                    test_node_cptr = 0
                    for test_node in (gen_paralog_node for gen_paralog_node in node_cladogram.children if gen_paralog_node.is_leaf() == False):
                        if test_node.state == "y":
                            test_node_cptr += 1
                    if test_node_cptr == len(node_cladogram.children):
                        node_cladogram.state = "y"  #ssi les noeuds internes fils ont tous y, alors le noeud ¨paralogue possède y
    
                elif node_cladogram.clade_type == "ortholog": #si orthologue
                    orth_leaf_gen = (child_orth_leaf for child_orth_leaf in node_cladogram.children if child_orth_leaf.is_leaf() == True) #feuilles directement branchées au noeud orthologue
                    orth_leaf_cptr = len([orth_leaf for orth_leaf in orth_leaf_gen if orth_leaf.state == "y"]) #nombre de feuilles ayant y
                    orth_paralog_gen = (node_child_orth for node_child_orth in node_cladogram.get_descendants() if node_child_orth.clade_type == "paralog") #on ne prend que les noeuds paralogues (paralog list est un générateur créé en intension)
                    test_n = len([paralog_node for paralog_node in orth_paralog_gen if paralog_node.state == "n"]) #nombre de noeuds paralogues sans y
                    if test_n == 0 and orth_leaf_cptr > 0: #si aucun des noeuds paralogues n'est "n" et si au minimum une des feuilles directement branchée au noeud orthologue a y
                        node_cladogram.state = "y" #alors le noeud orthologue à y
            
            for orth_node in (node for node in cladogram.traverse(strategy="postorder") if node.clade_type =="ortholog" and node.state == "y"): #remplit les lignées orthologues de y pour éviter d'avoir en synapomorphie deux noeuds inclus l'un dans l'autre
                for orth_node_child in orth_node.get_descendants():
                    orth_node_child.state = "y"
            
            #à la fin de l'assignation, les noeuds y sont récupérés dans synapomorphies_set
            for leaf_syn_y in (leaf_syn for leaf_syn in cladogram.iter_leaves() if leaf_syn.state == "y"): #générateur des feuilles y
                synapomorphies_set.add(find_synapomorphy(leaf_syn_y)) #liste des noeuds caractérisés par l'état (synapomorphie)
    
        #vérification du nombre de noeuds retenus:1 = synapomorphie (acceptation), 2 ou plus = homoplasie (rejet)
        if len(synapomorphies_set) > 1: #condition du rejet d'un état
            synapomorphy_test = str("rejected")
        results_test_dict[char_names][char_states_names] = [synapomorphy_test, synapomorphies_set]
            
        return synapomorphy_test, synapomorphies_set, results_test_dict, cladogram

    #chaque noeud de chaque cladogramme est numéroté et défini comme paralog/ortholog/apical
    for cladogram in cladogram_dict:
        cladogram, syn_dict = cladogram_label(cladogram, clade_number_option="yes", clade_type_option="yes") #label des noeuds: types et numéros
    
    #liste des caractères, de leurs états et leur contenu en taxons
    #taxa in state
    character_component_dict = {}
    for character, values in character_dict.items(): 
        character_states = {}
        character_states_count = 1
        for node in character.traverse(strategy="preorder"): #label des noeuds du cladogramme
            if node.is_leaf() == False  and node.is_root() == False:
                character_states["Character state #"+str(values)+"."+str(character_states_count)] = Tree.get_leaf_names(node) #dictionnaire de taxons porteurs de l'état dérivé
                character_states_count += 1
        character_component_dict["Character_"+str(values)] = character_states #liste des taxons dans le noeud
    #character cardinal
    character_cardinal_dict = {}
    for character, values in character_dict.items(): 
        character_states_count = 1
        for node in character.traverse(strategy="preorder"): #label des noeuds du cladogramme
            if node.is_leaf() == False  and node.is_root() == False:
                character_cardinal_dict["Character state #"+str(values)+"."+str(character_states_count)] = set(Tree.get_leaf_names(character)) #dictionnaire de taxons porteurs de l'état dérivé
                character_states_count += 1
    
    ##################################
    #TEST DES ETATS ET AJOUT DU STYLE#
    ##################################
                
    results_test_dict = {}
    for cladogram in cladogram_dict: #pour l'instant, un seul
        with open(prefix+".chartest", "a") as results_file: #hekate states tests
            results_file.write("#Character states tests")
            results_file.write("\n#Legend: #character.state / state accepted or rejected / node(s) characterized by the state (if one: synapomorphy)")
            results_file.write("\n")
    
        #pour chaque caractère
        for char_names, values in character_component_dict.items():
            results_test_dict[char_names] = {}
    
            #pour chaque état de caractère
            for char_states_names, taxa_in_state in values.items():           

                #test des états de caractères
                if character_cardinal_dict[char_states_names] == set(Tree.get_leaf_names(cladogram)): #s'il n'y a pas de données manquantes
                    synapomorphy_test, synapomorphies_set, results_test_dict, cladogram = synapomorphy_test_placement(cladogram, taxa_in_state, results_test_dict, char_names, char_states_names) #la fonction cherche quels sont les noeuds caractérisés par l'état dérivé pris en argument, ainsi que si l'état est accepté ou rejeté
    
                else: #s'il y a des données manquantes
                    pruned_cladogram = cladogram.copy(method='cpickle') #on fait une copie
                    pruned_cladogram.prune(list(character_cardinal_dict[char_states_names])) #sont élaguées les feuilles présentant des taxons portant des données manquantes (?)
                    pruned_cladogram, temp = cladogram_label(pruned_cladogram, clade_number_option="no", clade_type_option="yes") #label des noeuds: types et numéros
                    synapomorphy_test, synapomorphies_set, results_test_dict, pruned_cladogram = synapomorphy_test_placement(pruned_cladogram, taxa_in_state, results_test_dict, char_names, char_states_names) #la fonction cherche quels sont les noeuds caractérisés par l'état dérivé pris en argument, ainsi que si l'état est accepté ou rejeté
    
                    syn_list_extension = {} #construction du dictionnaires de synapomorphies définies par leurs feuilles (juste le nom)
                    for synapomorphy in synapomorphies_set:
                        if synapomorphy.is_leaf():
                            syn_list_extension[synapomorphy] = [synapomorphy.name]
                        elif synapomorphy.clade_type == "apical_node":
                            syn_list_extension[synapomorphy] = [leaf.name for leaf in synapomorphy.iter_leaves() if leaf.name in taxa_in_state] #les synapomorphies ne concernant que les noeuds apicaux de l'arbre élagué sont supprimés, et on n'en retir que les feuilles. Les noeuds apicaux seront vérifiés ensuite sur l'arbre complet
                        else:
                            syn_list_extension[synapomorphy] = [leaf.name for leaf in synapomorphy.iter_leaves()] #clé: synapomorphie; valeur: liste des feuilles
    
                    synapomorphies_set_2 = set() #set de feuilles/synapomorphies (avant recherche des feuilles liés à des noeuds apicaux)
                    synapomorphies_set = set() #set de synapomorphies réinitialisé
    
                    for syn_list in syn_list_extension.values(): #on dissocie les noeuds internes qui vont directement dans synapomorphies_set des feuilles qui vont dans synapomorphies_set_2 pour vérification
                        if len(syn_list) == 1: #forcément une feuille
                            synapomorphies_set_2.add(cladogram.get_leaves_by_name(syn_list[0])[0]) #isolement dans un autre set pour vérification
                        else:
                            synapomorphies_set.add(cladogram.get_common_ancestor([cladogram.get_leaves_by_name(leaf_name)[0] for leaf_name in syn_list])) #remplissage des noeuds internes par recherche de l'ancêtre commun de toutes les feuilles
    
                    for node in synapomorphies_set_2: #correction des feuilles y branchées aux noeuds apicaux 
                        if node.get_ancestors()[0].clade_type == "apical_node":
                            synapomorphies_set.add(node.get_ancestors()[0]) #feuille remplacée par son noeud ancêtre car il est apical
                        elif node.get_ancestors()[0].clade_type == "ortholog" and set([node_name.name for node_name in node.get_ancestors()[0].get_leaves()]) & set(taxa_in_state) == {node.name}: #si le noeud ancêtre ne contient que des taxons porteurs de données manquantes à part la feuille "node", alors assignation au noeud (ici, si l'intersection entre les taxons présents dans le caractère et les taxons présents dans le noeud ancêtre ne donnent que lo feuille)
                            synapomorphies_set.add(node.get_ancestors()[0]) #feuille remplacée par son noeud ancêtre car orthologue rempli de ?
                        else:
                            synapomorphies_set.add(node) #sinon ajout tel quel
    
                    results_test_dict[char_names][char_states_names] = [synapomorphy_test, synapomorphies_set] #on remet le synapomorphies_set mis à jour avec le cladogramme non-élagué
    
                #ajout des styles de feuilles / uniquement pour les pdf
                if pdf_files:
                    
                    synapomorphy_style = TreeStyle() #style de l'arbre
                    synapomorphy_style.show_leaf_name = True
                    synapomorphy_style.margin_top = 10
                    synapomorphy_style.margin_right = 10
                    synapomorphy_style.margin_left = 10
                    synapomorphy_style.margin_bottom = 10            
                    synapomorphy_style.show_scale = False
                    synapomorphy_style.show_leaf_name = False
                    synapomorphy_style.title.add_face(TextFace(char_states_names, fsize=3), column=0) #titre du pdf
        
                    for node1 in cladogram.iter_leaves():
                        if node1.name in taxa_in_state:
                            node1.add_feature("taxa_in_state", "in")
                        elif node1.name in character_cardinal_dict[char_states_names]:
                            node1.add_feature("taxa_in_state", "out")
                        else:
                            node1.add_feature("taxa_in_state", "?")
                            
                    synapomorphy_style.layout_fn = lstyle #ajout du layout modifiant les noeuds
                
                    #ajout des styles de noeuds internes
                    for node_style1 in cladogram.traverse(strategy="postorder"): 
                        if node_style1 in synapomorphies_set:
                            if synapomorphy_test == "accepted":
                                node_style1.set_style(nstyle_accept) #-ajout de l'information des feuilles portant l'état dérivé
                            else:
                                node_style1.set_style(nstyle_reject) #-ajout de l'information des feuilles portant l'état dérivé
                        else:
                            node_style1.set_style(nstyle_other)            
                    
                    #ajout des commentaires au pdf
                    synapomorphy_style.legend.add_face(TextFace("Synapomorphy test: "+synapomorphy_test, fsize=2), column=0)
                    synapomorphy_style.legend.add_face(TextFace("", fsize=2), column=0)
                    synapomorphy_style.legend.add_face(TextFace("Node(s): ", fsize=2), column=0)
                    node_line_pdf = str()
                    for node_text in sorted([syn_node.clade_label for syn_node in synapomorphies_set if not syn_node.is_leaf()]):
                        node_line_pdf += str(node_text)
                        node_line_pdf += " "
                    synapomorphy_style.legend.add_face(TextFace(node_line_pdf, fsize=2), column=0)
                    
                    node_line_pdf = str()
                    for node_text in sorted([syn_node.clade_label for syn_node in synapomorphies_set if syn_node.is_leaf()]):
                        node_line_pdf += str(node_text)
                        node_line_pdf += " "
                    synapomorphy_style.legend.add_face(TextFace(node_line_pdf, fsize=2), column=0)
                    
                    #enregistrement du pdf
                    cladogram.render(pdf_files+char_states_names+".pdf", tree_style=synapomorphy_style)
                    
                    
                #complétion du fichier texte
                with open(prefix+".chartest", "a") as results_file:
                    results_file.write("\n"+char_states_names+": "+synapomorphy_test)
                    syn_not_leaf = sorted([syn_node.clade_label for syn_node in synapomorphies_set if not syn_node.is_leaf()])
                    syn_leaf = sorted([syn_node.clade_label for syn_node in synapomorphies_set if syn_node.is_leaf()])
                    if len(syn_not_leaf) == 1:
                        results_file.write(", clade ")
                    elif len(syn_not_leaf) > 1:
                        results_file.write(", clades ")
                    for node_syn in syn_not_leaf:
                        results_file.write(str(node_syn)+", ") #str est placé ici pour permettre le tri des int dans syn_not_leaf
                    for node_syn in syn_leaf:
                        results_file.write(node_syn+", ")
                        
    #création du fichier avec les positions des synapomorphies
    for char_names, values1 in results_test_dict.items():
        for char_states_names, syn_test_set in values1.items(): #pour chaque couple test, node
            for node_set in syn_test_set[1]:
                if syn_test_set[0] == "accepted":
                    syn_dict[str(node_set.clade_label)]["accepted"].append([char_states_names.replace("Character state ","")])
                elif syn_test_set[0] == "rejected":
                    syn_dict[str(node_set.clade_label)]["rejected"].append([char_states_names.replace("Character state ","")])
    
    with open(prefix+".chartest_node", "w") as results_file_tree:
        results_file_tree.write("Synapomorphies by node:")
        results_file_tree.write("\n")
        for node_set in sorted((str(node_set) for node_set, values in syn_dict.items())):
            if syn_dict[node_set]["accepted"] or syn_dict[node_set]["rejected"]:
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

def child_triplet_generator(triplet1, triplet2): 
    """
        prend des triplet en arguments 
        Calcule les triplets secondaires générés par un couple de triplets
    """
    no_redundant_generated_triplet = set() #contient les triplet à renvoyer
    if len(set(triplet1.in_taxa | triplet1.out_taxa) & set(triplet2.in_taxa | triplet2.out_taxa)) == 2 and not triplet1.in_taxa == triplet2.in_taxa and triplet1.parent_triplets != triplet2.parent_triplets and triplet1 not in triplet2.parent_triplets and triplet2 not in triplet1.parent_triplets : #deux conditions pour entrer dans la boucle de génération des triplets: 1) avoir deux taxons en commun, 2) mais pas in-in, et 3) ne pas avoir les mêmes parents immédiats et 4) ne pas faire de combinaisons parent-enfant (correspond aux deux derniers and)
        
        if triplet1.out_taxa == triplet2.out_taxa and len(triplet1.in_taxa & triplet2.in_taxa) == 1: #génération du triplet in-ex
            generated_triplet = {triplet(triplet1.in_taxa ^ triplet2.in_taxa, triplet1.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level))}
        
        elif triplet1.out_taxa < triplet2.in_taxa and triplet2.out_taxa < triplet1.in_taxa: #génération des deux triplets paralogues
            generated_triplet = {triplet(triplet1.in_taxa, triplet2.in_taxa - triplet1.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level)), triplet(triplet2.in_taxa, triplet1.in_taxa - triplet2.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level))}
        
        elif len(triplet2.in_taxa & triplet1.out_taxa) == 1 and len(triplet1.in_taxa & triplet2.in_taxa) == 1: #génération des deux triplets orthologues
            generated_triplet = {triplet(triplet1.in_taxa, triplet2.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level)), triplet(triplet2.in_taxa ^ triplet1.in_taxa, triplet2.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level))}
        
        else: #génération des deux triplets orthologues (version 2)
            generated_triplet = {triplet(triplet2.in_taxa, triplet1.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level)), triplet(triplet2.in_taxa ^ triplet1.in_taxa, triplet1.out_taxa, parent_triplets = (triplet1, triplet2), primary_parent_triplets = triplet1.primary_parent_triplets | triplet2.primary_parent_triplets, gen_level = 1 + max(triplet1.gen_level, triplet2.gen_level))}
        
        for single_generated_triplet in generated_triplet:
            if not (single_generated_triplet.in_taxa, single_generated_triplet.out_taxa) in ((triplet_compare.in_taxa, triplet_compare.out_taxa) for triplet_compare in triplet1.primary_parent_triplets | triplet2.primary_parent_triplets): #test que chaque triplet généré soit bien inédit par rapport à la liste des triplets parents primaires
                no_redundant_generated_triplet.add(single_generated_triplet)
                print(single_generated_triplet) #soit zéro, soit un seul triplet, soit un tuple de deux triplet
        
    return no_redundant_generated_triplet #soit zéro, soit un seul triplet, soit un tuple de deux triplet

def describe_tree(dtree, outfile, nb=1, showtaxanames=False):
    """
    

    Parameters
    ----------
    dtree : ete3 tree
        DESCRIPTION.
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
    
    taxalist = taxa_list_extraction({dtree:0})
    repetitions = len(taxalist) - len(set(taxalist))
        
    for node in dtree.traverse("postorder"):
        if not node.is_leaf(): #si le noeud est un noeud interne
        
            totaln += 1
            i = len([child_node for child_node in node.get_children() if not child_node.is_leaf()]) #nombre de descendants noeuds internes
            j = len([child_node for child_node in node.get_children()]) #nombre de descendants #nombre de descendants noeuds
            
            if i >= 2: #paralog
                paralogn += 1
            elif i == 1: #ortholog
                orthologn += 1
            elif i == 0: #apical
                apicaln += 1
                
            if j == 2:#dichotomies vs polytomies
                dichotomies += 1
            else:
                polytomies += 1   
                
                if j in polytomiesd:#détail des polytomies
                    polytomiesd[j] += 1
                else:
                    polytomiesd[j] = 1
                
    with open(outfile, "a") as describefile:
        describefile.write("***************\n")
        describefile.write("Describe tree {}\n".format(str(nb)))
        #print(*terminals, sep = ", ") #ligne pour afficher les terminaux
        describefile.write("Number of nodes: "+str(totaln+cardinal)+"\n")
        describefile.write("    Terminals: "+str(cardinal)+"\n")
        describefile.write("    Internal nodes: "+str(totaln)+"\n")
        describefile.write("        Symmetric nodes: "+str(paralogn)+"\n")
        describefile.write("        Asymmetric nodes: "+str(orthologn)+"\n")
        describefile.write("        Apical nodes: "+str(apicaln)+"\n")
        describefile.write("\n")
        
        if cardinal > 2:
            describefile.write("Resolution: "+str(round(((totaln-1)/(cardinal-2))*100,5))+"%\n") #on ne compte pas la racine
        else:
            describefile.write("Resolution: NA\n") #on ne compte pas la racine
        
        describefile.write("Dichotomies: {}; polytomies: {}\n".format(str(dichotomies),str(polytomies)))
        describefile.write("Repetitions:{}\n".format(str(repetitions)))
        
        if polytomies > 0:
            describefile.write("Polytomies - details (level: number of occurences)\n")
            for poly in sorted(polytomiesd.keys()):
                describefile.write("{} : {}\n".format(str(poly), str(polytomiesd[poly])))
       
        if showtaxanames:
            describefile.write("\nTaxa list:\n")
            for taxa in taxalist:
                describefile.write(taxa+"\n")
        
        
def describe_forest(character_dict, prefix, showtaxanames=False):

    i = 0
    for dtree in character_dict.keys():
        i += 1
        describe_tree(dtree, prefix+".dt", nb=i, showtaxanames=showtaxanames)


