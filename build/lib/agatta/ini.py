# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:46:58 2020

@author: Valentin Rineau
"""

from re import findall, search, compile
from tkinter import filedialog
from tkinter import Tk
from os import path
from collections import defaultdict
import warnings
import csv

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=SyntaxWarning)
    from ete3 import Tree

def helper(option):
    """
    

    Parameters
    ----------
    option : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if option == "analysis":
        print("Help block line for analysis command")
    elif option == "tripdec":
        print("Help block line for tripdec command")
    elif option == "ri":
        print("Help block line for ri command")
    elif option == "chartest":
        print("Help block line for chartest command")
    elif option == "convert":
        print("Help block line for convert command")
    elif option == "fp":
        print("Help block line for fp command")
    elif option == "consensus":
        print("Help block line for consensus command")
    elif option == "describetree":
        print("Help block line for describetree command")
    elif option == "standardisation":
        print("Help block line for standardisation command")
    elif option == "hmatrix":
        print("Help block line for hmatrix command")

def hmatrix(nexus_file_name):
    """
        prend en argument un (nexus) path
        
        fonction prenant un fichier nexus PAUP en entrée
        extraction de l'outgroup et de la matrice taxons/caractères
        construction d'arbres de caractères /!\ attention, ici tout
        est ordonné linéairement, par ex: 0->1->2->3 ou 0<-1<-2->3 /!\
        polarisation par extra-groupe à spécifier dans le nexus
    """
    
    #déclaration des variables
    matrix_lines = {}
    taxa_list = []
    character_dict = {}
    char_number = 0
    
    #lecture du fichier nexus
    lecture_var = ""
    
    with open(nexus_file_name, "r") as nexus_file:
        for line in nexus_file:
            
            if search("taxlabels", line): #activation du mode récupération des taxons à partir de la ligne suivante
                lecture_var = "taxa"
    
            elif search("matrix", line): #activation du mode récupération de la matrice à partir de la ligne suivante
                lecture_var = "matrix"
                
            elif lecture_var == "taxa" and search("(\t)((.)*)", line): #mode récupération des taxons
                taxa_list += search("(\t)((.)*)", line).group(2).split(" ")
                if search(";", taxa_list[-1]):
                    taxa_list[-1] = taxa_list[-1][:-1] #supression du point-virgule
                    lecture_var = ""
    
            elif lecture_var == "matrix" and search("(\t)((.)*)( )((\d)+)", line): #mode récupération de la matrice
                matrix_lines[search("(\t)((.)*)( )((\d)+)", line).group(2)] = search("(\t)((.)*)( )((\d)+)", line).group(5)
                if search(";", line):
                    lecture_var = ""
                    
            elif search("(dimensions nchar=)((\d)+)", line): #nombre de caractères
                char_number = int(search("(dimensions nchar=)((\d)+)", line).group(2))
                                                                           
            elif search("(Outgroup )((.)*)(;)", line):
                outgroup_name = search("(Outgroup )((.)*)(;)", line).group(2) #nom de l'extra-groupe, qui permet 
                                                                            #la polarisation de chacun des caractères
            else:
                continue
    
    char_count = 0
    
    while char_count != char_number:
        
        t = Tree() #création de l'arbre
        
        #enracinement des feuilles à la racine
        state_number = int(matrix_lines[outgroup_name][char_count])
        state_number_count = int(matrix_lines[outgroup_name][char_count])
        
        for taxa in taxa_list: #branchement des taxons à la racine
            if int(matrix_lines[taxa][char_count]) == state_number:
                t.add_child(name=taxa)
    
        #ajout d'un noeud interne branché à la racine pour l'état outgroup - 1
        if state_number - 1 >= 0:
            ancestor_n = t.add_child(name = "ancestor_"+str(state_number - 1)) # génération du noeud interne
            
            #ajout des noeuds internes orthologues suivants (state_number_count décroissant)
            while state_number_count >= 0: #pour chaque état jusqu'à 11 (à changer), ajout d'un noeud interne (état de caractère)
                state_number_count -= 1
                for taxa in taxa_list: #ajout des feuilles
                    if int(matrix_lines[taxa][char_count]) == state_number_count: 
                        ancestor_n.add_child(name=taxa) 
        
                ancestor_n = ancestor_n.add_child(name="ancestor_"+str(state_number)) # Adds a new child to the current tree root
            
        state_number_count = state_number
        
        #ajout d'un noeud interne branché à la racine pour l'état outgroup + 1
        if state_number + 1 <= 10:
            ancestor_n = t.add_child(name = "ancestor_"+str(state_number + 1)) # génération du noeud interne de l'état outgroup + 1
        
            #ajout des noeuds internes orthologues suivants (state_number_count croissant)
            while state_number_count <= 10: #pour chaque état jusqu'à 11 (à changer), ajout d'un noeud interne (état de caractère)
                state_number_count += 1
                for taxa in taxa_list: #ajout des feuilles
                    if int(matrix_lines[taxa][char_count]) == state_number_count: 
                        ancestor_n.add_child(name=taxa) 
        
                ancestor_n = ancestor_n.add_child(name="ancestor_"+str(state_number)) # Adds a new child to the current tree root
                        
            char_count += 1 #boucle while continue jusqu'à 100
            character_dict[t] = str(char_count)
    
    for character, value in character_dict.items():
        character.prune(taxa_list)
        
    return character_dict #dictionnaire de caractères hiérarchiques déduits de la matrice du nexus

def character_extraction(infile, taxa_replacement=False): 
    """
        Fonction demandant un fichier et renvoyant un dictionnaire de caractères 
        (si taxa_replacement="yes", les symboles sont remplacés par les 
        taxons complets)
        /!\ NE MARCHE QU'AVEC LES FICHIERS .3IA + FORMAT NEWICK /!\
    """
    
    print("Loading character trees")
    
    character_dict = {} #dictionnaire ou sont enregistrés l'ensemble des arbres de caractères du fichier "pathmessage"
    taxa_dict = {} #dictionnaire de conversion des taxons/symboles

    if infile == None:
        root = Tk()
        infile =  filedialog.askopenfilename(title = "Select file containing newick character trees",filetypes = (("all files","*.*"),("tree files","*.tre"),("3ia files","*.3ia"),("nexus files","*.nex")))
        root.withdraw()

    fileName, fileExtension = path.splitext(infile)

    with open(infile, "r") as file_tree :
        a=1
        line_nb = 0
        for line in file_tree: 
            line_nb += 1
            for character_newick in findall(r"\([^ \t\n\r\f\v]+;", line):
                if character_newick:
                    try:
                        character_dict[Tree(character_newick)] = a #dictionnaire d'arbres (indicé avec la variable a)
                        a += 1
                    except:
                        print("Line {}: Broken newick structure.".format(str(line_nb)))                    
                    
            if taxa_replacement and fileExtension == ".3ia" and search(r"\s=\s", line): #recherche du taxa bloc en 3ia en remplacement des noms si taxa_replacement="yes"
                taxa_dict[line.split(" = ")[0].split()[0]] =  line.split(" = ")[1].strip()
                
    if taxa_dict:
        for cladogram in character_dict.keys():
            for leaf in cladogram.iter_leaves():
                if taxa_dict[leaf.name]:
                    leaf.name = taxa_dict[leaf.name]
                    
    print("{} characters loaded".format(str(len(character_dict))))
    
    return character_dict


def taxa_list_extraction(character_dict): 
    """
        extrait un ensemble de taxons à partir d'un dictionnaire 
        d'arbres de caractères
    """
    taxa_list = set() #ensemble des taxons utilisés
    for character in character_dict.keys(): #génération de la matrice PAUP avec taxons
        taxa_list.update(set(Tree.get_leaf_names(character)))
    return taxa_list

def taxa_triplet_extraction(triplet_dict):
    """
    prend un dictionnaire de triplets et renvoie une liste de taxons
    """
    
    taxa_list = set() #ensemble des taxons utilisés
    
    for triplet in triplet_dict.keys(): #génération de la matrice PAUP avec taxons
        taxa_list.update(triplet.out_taxa)
        taxa_list.update(triplet.in_taxa)
        
    return taxa_list
    

def taxa_to_numbers(character_dict): 
    """
        prend un ensemble de taxons et renvoie les taxons avec un identifiant 
        (nombre associé)
    """
    
    taxa_list = taxa_list_extraction(character_dict)
    taxa_dict = {}
    code_dict = {}
    taxa_convers = []
    taxa_cptr = 1
    for taxa1 in taxa_list: #pour chaque taxon
        taxa_dict[taxa1] = taxa_cptr #integer labels
        code_dict[taxa_cptr] = taxa1 #integer labels
        taxa_convers.append(str(taxa1)+" "+str(taxa_cptr))
        taxa_cptr += 1
        
    return taxa_dict, code_dict, taxa_convers #1) dictionnaire avec taxons et leurs numéros associés; 2) même chose en bloc texte, sous forme de liste

def taxa_to_lisbeth_numbers(taxa_list): 
    """
        prend un ensemble de taxons et renvoie les taxons avec un identifiant 
        (nombre associé)
    """
    taxa_dict = {}
    taxa_txt_bloc = []
    taxa_cptr1 = 65
    taxa_cptr2 = 65
    for taxa1 in sorted(taxa_list): #pour chaque taxon
        taxa_dict[taxa1] = chr(taxa_cptr1) + chr(taxa_cptr2)
        taxa_txt_bloc.append(" "+chr(taxa_cptr1) + chr(taxa_cptr2)+" = "+taxa1)
        taxa_cptr2 += 1
        if taxa_cptr2 == 91:
            taxa_cptr1 += 1 #la première lettre est incrémentée
            taxa_cptr2 = 65 #la deuxième lettre redevient A     
        
    return taxa_dict, taxa_txt_bloc #1) dictionnaire avec taxons et leurs numéros associés; 2) même chose en bloc texte, sous forme de liste

def phylo_to_areagram(tree_file, biogeo_tab, prefix="agatta_standardisation", 
                      verbose=False):
    def biogeo_table(biogeo_tab, verbose=False):
        """
            fonction prenant un fichier contenant des arbres newick ou un chardict
            table de conversion entre taxons et aires biogéographiques
        """
        
        #déclaration des variables
        biogeo_dict = {}
        taxa = set() #set total des taxons
        areas = set() #set total des aires
        MAST = set() #set des clefs correspondant aux MAST
        taxrep = set() #set des items correspondant à de la répétition de terminaux
    
        
        with open(biogeo_tab, "r") as bt_file:
            for line in bt_file:
                regex = compile(r'[\n\r\t]')
                line = regex.sub("",line)
                linesplt = line.split(":")
                
                if linesplt[0] in taxa: #si taxon déjà présent => MAST
                    biogeo_dict[linesplt[0]].append(linesplt[1])
                    MAST.add(linesplt[0])
                else:
                    biogeo_dict[linesplt[0]] = [linesplt[1]] #clef: taxon, item, aire
                
                if linesplt[1] in areas: #si aire déjà présente => répétition
                    taxrep.add(linesplt[1])
                    
                taxa.add(linesplt[0])
                areas.add(linesplt[1])
                
        if verbose:
            print("{} terminal taxa, {} terminal areas".format(str(len(taxa)),str(len(areas))))
            print("{} repeated taxa, {} MAST detected".format(str(len(taxrep)),str(len(MAST))))
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
    
    if nb_char > 1 :
        print("Standardising {} characters".format(str(nb_char)))
    else:
        print("Standardising 1 character")
    
    #pour chaque arbre
    with open(prefix+".std", "w") as area_file:
        for phylogeny, index in character_dict.items():
            
            phylogeny2 = phylogeny.copy() #copie de l'arbre
            
            #supression des MAST
            for leaf in phylogeny.get_leaves():
                if leaf.name in MAST: #si MAST on supprime la feuille
                    leaveswmast = list(set(phylogeny2.get_leaf_names()) - MAST)
                    phylogeny2.prune(leaveswmast)
            
            #pour chaque terminal taxon, remplacement par aire * sauf MAST
            for leaf in phylogeny.get_leaves():
                if not leaf.name in MAST: #si MAST on supprime la feuille
                    for leafnode in phylogeny2.get_leaves_by_name(leaf.name):
                        leafnode.name = biogeo_dict[leaf.name][0]
                                    
            #fin, ajouter le nouvel arbre à la liste
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
        
def matrix_to_trees(path, prefix=False, chardec=False, verbose=False):
    """
    
    Parameters
    ----------
    path : path csv file string
    verbose : if verbose, show managed duplications

    Returns 
    -------
    character_dict : TYPE
        DESCRIPTION.
    polytomies : dictionary
        Return a dictionary character number: list of repeated taxa names.

    """
    
    print("Loading hierarchical matrix")
    
    character_dict = dict() #arbres sans polytomies
    temp_character_dict = dict() #arbres avec polytomies (raw data)
    error_dict = defaultdict(list)
    
    #lecture de la matrice
    with open(path,'rt')as f:
      data = csv.reader(f,delimiter=';')
      hmatrix = list(data)
      
    print("Hierarchical matrix loaded")
    print("Treefication of the hierarchical matrix.")
    
    #construction des caractères non remplis
    i = 1     
    treeliststr = hmatrix[0]
    del treeliststr[0]
       
    for char in treeliststr:
        temp_character_dict[Tree(char+";")] = str(i)
        i += 1
        
    taxalist = [hmatrix[ncol][0] for ncol in range(1,len(hmatrix))]
    polytomies = defaultdict(list)
    
    #option pour décomposer en composantes
    if chardec:
        temp_character_dict_binary = dict()
        
        for char, ind in temp_character_dict.items(): #pour chaque caractere
            
            #pour chaque état de caractere non racine
            for state in char.traverse(strategy="levelorder"):
                if state.is_leaf() and not state.get_ancestors()[0].is_root():
                    
                    statetree = char.copy()
                    
                    for isolstate in statetree.traverse(strategy="levelorder"):
                        if not isolstate.is_leaf() and not isolstate.is_root():
                            if not isolstate.get_children()[0].name == state.name:
                                isolstate.delete()
                    
                    temp_character_dict_binary[statetree] = str(ind)+"."+str(state.name)
                    
        temp_character_dict = temp_character_dict_binary
                
        
    #remplissage des caractères
    for char, ind in temp_character_dict.items(): #pour chaque caractere
    
        #marquage des branches à supprimer à la fin    
        delnodes = [] #liste des feuilles à supprimer
        for leaf in char.iter_leaves():        
            delnodes.append(leaf)
            leaf.name = "_agatta_charstate_"+leaf.name #évite le bug ou même nom pour état et pour taxon
    
        #pour chaque taxon
        for taxa in taxalist:
            
            if chardec:
                ind2 = ind.split(".")[0]
            else:
                ind2 = ind
            
            charstates = hmatrix[taxalist.index(taxa)+1][int(ind2)].split(",") #état ou placer le taxon
            
            if charstates[0] != "?": #si pas donnée manquante
            
                # #construction d'un dictionnaire numéro de caractere vs taxons répétés
                # if len(charstates) > 1:
                #     polytomies[ind2].append((taxa,len(charstates)))
                
                #branchement
                for charstate in charstates:
                    
                    #si la matrice est mal faite, une erreur est capturée ici, sinon, branchement du taxon
                    try:
                        branchnode = char.get_leaves_by_name("_agatta_charstate_"+charstate)[0].get_ancestors()[0] #noeud ou brancher le taxon
                        branchnode.add_child(name=taxa)
                        
                    except:
                        error_dict[ind2].append(taxa) #indice du caractere / taxon
                        
        #supression des branches d'états
        for delnode in delnodes:
            delnode.delete()    
        
        
    #gestion des répétitions de taxons (polymorphisme)
    #polytomchar = [charnumber for charnumber, taxalist in polytomies.items()]
        
    for char, value in temp_character_dict.items():
        
        # #si polytomie détectée
        # if value in polytomchar and method != "No":
            
        #     if verbose:
                
        #         print("Polytomies detected for character number {}. Automatic polytomies deletion performed.".format(str(value)))
            
        #     treelist = del_replications(char, method, False)
            
        #     for treepoly in treelist:
        #         treepoly.ladderize()
        #         character_dict[treepoly] = str(value)
    
        # #si pas de polytomies
        # else:
            
        char.ladderize()
        character_dict[char] = str(value)
           
    #mode verbose
    if verbose:
        if error_dict:
            
            print("Error in the input matrix\n")
            
            for charnum, taxalist in error_dict.items():
                print("Character number {}:\n".format(str(charnum)))
                for taxa in taxalist:
                    print(" - "+taxa)
                
    #save file
    if prefix:
        # with open(prefix+".tre", "w") as treefile:
        #     for char in character_dict.keys():
        #         treefile.write(char.write(format=9)+"\n")

        with open(prefix+".hmatrix", "w") as logfile:
            for char, charnum in character_dict.items():
                logfile.write(str(charnum)+" : "+char.write(format=9)+"\n")
            
            # logfile.write("[POLYTOMIES]\n")
            # for char, polytuple in polytomies.items():
            #     logfile.write("Character number {}:\n".format(char))
            #     for taxaname, nrepet in polytuple:
            #         logfile.write("   Taxa {} repeated {} times\n".format(taxaname,nrepet))

    print("{} characters computed from the matrix".format(str(len(character_dict))))        

    return character_dict, polytomies, error_dict #polytomies: dictionnaire numéro de caractere : [(taxaname,nombre de répétitions),(),()]

        
        

        
        

