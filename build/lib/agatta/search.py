# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:27 2020

@author: Valentin Rineau
"""

import os, platform
from random import choice
from .ini import taxa_list_extraction
import warnings
import time

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
    from ete3 import Tree


def search_pipeline(path_infile, software_path, software="paup"):
    """
    

    Parameters
    ----------
    path_infile : string
        path of the nexus or tnt file.
    software_path : string
        path of the software chosen, paup or tnt.
    software : string, optional
        paup or tnt. The default is "paup".

    Raises
    ------
    Exception
        exception if os is windows.

    Returns
    -------
    nothing. run the analysis with specified file and software

    """
    
    print("Running analysis on "+software+" software.")
    start = time.time()
    
    ostype = platform.system() #détection de l'os
        
    if ostype == "Windows":
        if software == "paup":
            os.system("paup "+path_infile)

        elif software == "tnt":
            raise Exception("Windows is currently unsupported by Agatta. Currently supported systems: Linux, Mac")
        
    else: #si l'os est linux ou mac
        if software == "paup":
            os.system(software_path+" -n "+path_infile+" > /dev/null")
            
        elif software == "tnt":
            os.system(software_path+" proc "+path_infile)
            
    print("Three-item analysis done")

    end = time.time()
    time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print("elapsed time (3ia analysis): {}".format(time_cptr))
        

def rear_taxa(tree1): 
    """
        retourne le même arbre dont la position d'un taxon à changé
        prend en argument un seul arbre (SPR à une feuille)
    """
    character_dict = {}
    character_dict[tree1.copy(method="newick")] = str(1)
    tree_copy_test = tree1.copy(method="newick")
    rf = 0 #score de similarité entre l'arbre de départ et l'arbre final
    
    while rf == 0: #tant que les arbres (newick strings) sont identiques
        tree2 = tree_copy_test.copy(method="newick") #réinitialisation de l'arbre
        rand_node = Tree() #la position à laquelle sera placée le taxon à déplacer (en groupe frère du noeud rand_node)
        rand_taxa = "" #le taxon à déplacer
        
        while rand_node.name == rand_taxa: #vérifie que le taxon de départ et la position finale ne sont pas les mêmes
            rand_node = choice([node for node in tree2.traverse(strategy="postorder")]) #noeud pris au hasard (position finale)
            rand_taxa = choice([taxa for taxa in taxa_list_extraction(character_dict)]) #taxon pris au hasard 
                                                                                        #(le taxon qui va être déplacé)
        regraft_node = rand_node.copy(method="newick")
        tree2.search_nodes(name=rand_taxa)[0].get_ancestors()[0].delete()
        tree2.search_nodes(name=rand_taxa)[0].delete() #on supprime le taxon à changer de place
        
        if regraft_node.search_nodes(name=rand_taxa):
            regraft_node.search_nodes(name=rand_taxa)[0].get_ancestors()[0].delete()
            regraft_node.search_nodes(name=rand_taxa)[0].delete() #on supprime le taxon à changer de place
        
        if len([leaf for leaf in rand_node.iter_leaves()]) == 1:
            leaf_name = rand_node.get_leaves()[0].name
            rand_node.get_leaves()[0].name = "new_ancestor"
            new_ancestor = rand_node.get_leaves()[0]
            new_ancestor.add_child(name=leaf_name)
            new_ancestor.add_child(name=rand_taxa)
        
        else:
            for node in rand_node.get_children(): #supression du noeud
                node.detach()
            new_ancestor = rand_node.add_child(name="new_ancestor") #regraft node + leaf
            new_ancestor.add_child(regraft_node)
            rand_node.add_child(name=rand_taxa)
        
        for node in tree2.traverse(strategy="levelorder"):
            if len(node.get_children()) == 1:
                node.children[0].delete()
        rf, max_rf, common_leaves, parts_t1, parts_t2, set1, set2 = tree2.robinson_foulds(tree_copy_test)
    
    return tree2, tree2.search_nodes(name=rand_taxa)[0] #retourne le même arbre dont la position d'un taxon à changé


def tripletscore(triplet_dict, stree):
    """
    Calcule le score d'un arbre dichotomique pour analyse bandb
    """
    score = 0
    
    for triplet, FW in triplet_dict.items():
                
        tripin = list(triplet.in_taxa)
        a = stree.get_leaves_by_name(tripin[0])
        b = stree.get_leaves_by_name(tripin[1])
        c = stree.get_leaves_by_name(list(triplet.out_taxa)[0])
        
        
        if a and b and c: #si les trois feuilles existent dans l'arbre
        
            common =  stree.get_common_ancestor([a[0],b[0]])
            cc = common.get_leaves_by_name(list(triplet.out_taxa)[0])
                            
            if cc: #si fail
                score += FW*2
            else: #si triplet dans l'arbre
                score += FW
                    
    return score


def bandb(leaves, triplet_dict, base_tree=False, optimal_score=False, optimal_tree_list=[]):
    """
    

    Parameters
    ----------
    leaves : TYPE
        liste de feuilles.
    triplet_dict : TYPE
        DESCRIPTION.
    base_tree : TYPE, optional
        DESCRIPTION. The default is False.
    optimal_score : TYPE, optional
        DESCRIPTION. The default is False.
    optimal_tree_list : TYPE, optional
        DESCRIPTION. The default is [].

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    l = leaves.copy()
    
    if not optimal_score:
        optimal_tree = Tree()
        optimal_tree.populate(len(leaves), leaves)
        optimal_tree_list=[optimal_tree]
        optimal_score = tripletscore(triplet_dict, optimal_tree)
    
        base_tree = Tree() #construction du 2is de départ
        base_tree.add_child(name=l[0])
        base_tree.add_child(name=l[1])
        
        l.pop(0)
        l.pop(0)
        
        #print("debut")
        #print(float(optimal_score))
        #print(optimal_tree)
        #print(base_tree)
                    
    leaf = l[0]
    l.pop(0)
    
    for node in base_tree.traverse(): #pour chaque noeud, on créé un noeud d'insertion dans la branche qui mène vers lui
        
        #print("loop")
        #étape 1 - réarrangements
        if node.is_root(): #cas spécial de la nouvelle feuille à brancher en groupe frère du reste

            newtree = Tree()
            insert = base_tree.copy()
            newtree.add_child(insert)
            newtree.add_child(name=leaf)
            
            #print("root")
            
        else: #autres placements
            
            newtree = base_tree.copy() #arbre local
                        
            if node.is_leaf(): #chercher le noeud ancêtre
                newnode = newtree.get_leaves_by_name(node.name)[0]
            else:
                newnode = newtree.get_common_ancestor([l.name for l in node.get_leaves()])
                
            anc = newnode.get_ancestors()[0] #pour chaque noeud, prendre l'ancêtre
            insertnode = anc.add_child(name="internal") #y brancher un noeud
            insert = newnode.copy()
            newnode.detach() #supprimer le sous arbre
            insertnode.add_child(insert) #brancher le sous-arbre inclus dans ancêtre sur la feuille
            insertnode.add_child(name=leaf) #ajouter la nouvelle feuille

        #étape 2 - calcul du score et récursivité
        tripscore = tripletscore(triplet_dict, newtree) #calcul du score de l'arbre
            
        #print(str(float(tripscore))+" - "+str(float(optimal_score)))
        #print(newtree)
    
    
        if tripscore > optimal_score: #si l'arbre à déjà un score sup, on arrête
            
            #print("arrêt recursion")
            
            return [optimal_score, optimal_tree_list] #fail, stop this search
                            
        #print(l)
        
        if l: #si la récursion n'est pas finie
            bandb_result = bandb(l, triplet_dict, newtree, optimal_score, optimal_tree_list)
                            
            if bandb_result[0] < optimal_score: #si la récursion donne un résultat
                optimal_score=bandb_result[0] #nouveau score optimal
                optimal_tree_list=bandb_result[1] #nouvel arbre optimal  
                
            if bandb_result[0] == optimal_score:
                
                add_tree = True
                
                for rftree in optimal_tree_list:
                    for rftree2 in bandb_result[1]:
                        if rftree.robinson_foulds(rftree2)[0] == 0:
                            add_tree = False
                
                if add_tree:
                    optimal_tree_list.append(bandb_result[1])
                    
        else: #si la récursion est finie et que le résultat est meilleur
            #print("fin récursion")
            
            if tripscore == optimal_score:
                
                add_tree = True
                
                for rftree in optimal_tree_list:
                    if rftree.robinson_foulds(newtree)[0] == 0:
                        add_tree = False
                
                if add_tree:
                    optimal_tree_list.append(newtree)
                    
            if tripscore < optimal_score: #si la récursion donne un résultat
                optimal_score=tripscore #nouveau score optimal
                optimal_tree_list=[newtree] #nouvel arbre optimal  
                
            if tripscore == optimal_score:
                
                add_tree = True
                
                for rftree in optimal_tree_list:
                    if rftree.robinson_foulds(newtree)[0] == 0:
                        add_tree = False
                
                if add_tree:
                    optimal_tree_list.append(newtree)
                    

    return optimal_score, optimal_tree_list
    