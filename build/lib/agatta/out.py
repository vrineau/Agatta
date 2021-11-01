# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 10:47:17 2020

@author: Valentin Rineau
"""

from tkinter import filedialog
from tkinter import Tk
from .ini import character_extraction, phylo_to_areagram, hmatrix
from .ini import taxa_list_extraction, taxa_to_numbers, taxa_triplet_extraction
from .analysis import triplet_extraction
from .analysis import del_replications_forest
from .analysis import main_tripdec
from .analysis import rep_detector
from .interpret import RI, character_states_test
from .interpret import constrict, rcc
from .search import search_pipeline
from .search import bandb
import datetime
import time
import os
from ete3 import Tree


def triplet_nexus_file(tripletdict, character_dict, weighting, analysis, 
                       prefix=None, nrep=1000,logfile=False): 
    """
        Produit un fichier nexus à partir d'un dictionnaire triplets/poids
        Prend en argument le dictionnaire de caractères/triplets de triplet_decomposition
    """

    #écriture des poids
    count_zeroweights = 0 #ligne contenant le nombre de triplets à poids nul (arrondi)

    if weighting != "NW": #lignes des poids sauf pour NW
        i = 1 #triplet identifier    
        triplet_list_FW = str() #construction de la liste de pondérations
        deltriplets = dict()
    
        for trip, FW in tripletdict.items():    

            #management des pondérations en utilisant la pondération maximale de PAUP (automatique)
            if weighting == "FW" or weighting == "FWNL" or weighting == "MW": #tous ceux-ci sont des fractions
                if character_dict:
                    FWdec = len(str(len(character_dict)))
                else: FWdec = 2
                nexw = round((float(FW)),8-FWdec)
                if nexw == 0:
                    deltriplets[trip] = 0  #dans le cas ou l'arrondi donne zéro
                    count_zeroweights += 1
                elif nexw.is_integer():
                    triplet_list_FW += str(int(FW)) #sinon on insère l'arrondi (1.0 => 1)
                else:
                    triplet_list_FW += str(nexw) #sinon on insère l'arrondi
            else:
                nexw = 1
                triplet_list_FW += str(int(FW)) #NW, AW, UW (entiers)
    
            if nexw != 0:
                triplet_list_FW += ":"
                triplet_list_FW += str(i)
                triplet_list_FW += ", "
                i += 1

        for trip in deltriplets.keys():
            tripletdict.pop(trip)

    if character_dict:
        taxa_list = taxa_list_extraction(character_dict)
    else:
        taxa_list = taxa_triplet_extraction(tripletdict)
    
    #génération du fichier nexus
    if prefix==None:
        root = Tk()
        prefix =  filedialog.asksaveasfilename(title = "Save three-taxon analysis PAUP nexus file", filetypes = (("nexus files","*.nex"),("all files","*.*")))
        root.withdraw()

    with open(prefix+".nex", "w") as nexus_file:
        nexus_file.write("#NEXUS")
        nexus_file.write("\nbegin data;")
        nexus_file.write("\nDimensions ntax={} nchar={};".format(str(1+len(taxa_list)),str(len(tripletdict))))
        nexus_file.write("\nFormat symbols=\"0 1\" missing=?;")
        nexus_file.write("\n")
        nexus_file.write("\nMatrix")
        nexus_file.write("\n")

        #création de la matrice
        for taxa_name in taxa_list: # remplissage de la matrice PAUP avec les triplets
            nexus_file.write("\n"+taxa_name+" ")
            newline = ""
        
            #for each triplet
            for triplet_column in tripletdict.keys():
                if taxa_name not in triplet_column.out_taxa | triplet_column.in_taxa:
                    newline += "?"
                elif taxa_name in triplet_column.in_taxa:
                    newline += "1"    
                else:
                    newline += "0"
           
            nexus_file.write(newline) 
           
        nexus_file.write("\nroot "+"0"*len(tripletdict)) #all 0 root
            
        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")
        nexus_file.write("\n")
        nexus_file.write("\nBegin Paup;")
        
        #écriture des poids
        if weighting != "NW": #lignes des poids sauf pour NW
            nexus_file.write("\nwts")
            nexus_file.write("\n"+triplet_list_FW[:-2])     
        
        nexus_file.write("\n;")
        nexus_file.write("\noutgroup root /only;")
        
        if logfile:
            nexus_file.write("\nlog file={}_paup.log;".format(prefix))
        
        nexus_file.write("\nset maxtrees = 1000 increase=auto;")
        
        if analysis == "bandb":
            nexus_file.write("\nbandb;")
        elif analysis == "heuristic":
            nexus_file.write("\nhsearch addseq=random nreps={};".format(nrep))
            
        nexus_file.write("\nroottrees;")
        nexus_file.write("\nsavetrees /file={}.tre format=newick;".format(prefix))
#        nexus_file.write("\ncontree all / majrule=no strict=yes treefile={}_constrict.tre format=newick;".format(prefix))

        if logfile:
            nexus_file.write("\nlog stop;") #end of log save

        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")
        
        if count_zeroweights > 0:
            print("Warning, {} zero weight triplets deleted due to rounding".format(count_zeroweights))


def triplet_tnt_file(tripletdict, character_dict, weighting, analysis, prefix=None, nrep=1000, 
                       logfile=False): 
    """
        Produit un fichier pour tnt à partir d'une liste de triplets pondérés
        Prend en argument le dictionnaire de caractères/triplets de triplet_decomposition
    """    
    count_zeroweights = 0 #ligne contenant le nombre de triplets à poids nul (arrondi)            
        
    if weighting != "NW": #lignes des poids sauf pour NW
        triplet_list_FW = "ccode "
        wfactor = 1000/float(max(tripletdict.values()))
        itrip = 1
        deltriplets = dict()
        
        #construction de la matrice
        #if rescale weight needed
        if weighting == "FW" or weighting == "FWNL" or weighting == "MW" or max(tripletdict.values()) > 1000: #tous ceux-ci sont des fractions
            for triplet_column, FW in tripletdict.items(): #for each triplet
    
                nexw = int(float(FW)*wfactor)
                if nexw == 0:
                    deltriplets[triplet_column] = 0  #dans le cas ou l'arrondi donne zéro
                    count_zeroweights += 1
                else:
                    triplet_list_FW += " /{} ".format(str(nexw)) #sinon on insère l'arrondi
                    triplet_list_FW += str(itrip)
                    
                    itrip += 1
        
        #if rescale weight not needed
        else:
            for triplet_column, FW in tripletdict.items(): #for each triplet

                nexw = 1
                triplet_list_FW += " /{} ".format(str(FW))
                triplet_list_FW += str(itrip)
                
                itrip += 1

        triplet_list_FW += "*;" #fin de la ligne ccode  
        for trip in deltriplets.keys():
            tripletdict.pop(trip)

    #génération du fichier nexus
    if prefix==None:
        root = Tk()
        prefix =  filedialog.asksaveasfilename(title = "Save three-taxon analysis PAUP nexus file", filetypes = (("nexus files","*.nex"),("all files","*.*")))
        root.withdraw()

    if character_dict:
        taxa_list = taxa_list_extraction(character_dict)
    else:
        taxa_list = taxa_triplet_extraction(tripletdict)
    
    #construction de la chaine de caractère a écrire
    tntstring = "xread\n{} {}".format(str(len(tripletdict)),(str(1+len(taxa_list))))
    tntstring += "\n"
        
    #création de la matrice
    for taxa_name in taxa_list: # remplissage de la matrice PAUP avec les triplets
        tntstring += "\n"+taxa_name+" "
        newline = ""
    
        #for each triplet
        for triplet_column in tripletdict.keys():
            
            if taxa_name not in triplet_column.out_taxa | triplet_column.in_taxa:
                newline += "?"
            elif taxa_name in triplet_column.in_taxa:
                newline += "1"    
            else:
                newline += "0"
       
        tntstring += newline
        
    tntstring += "\nroot "+"0"*len(tripletdict) #all zeros root line
    
    tntstring += "\n;"
    tntstring += "\n"
    
    #construction de la ligne des poids
    if weighting != "NW": #lignes des poids sauf pour NW
        tntstring += "\n"+triplet_list_FW

    tntstring += "\n;"
    
    tntstring += "\n"
    
    if logfile:
        tntstring += "\nlog {}.tnt.log;".format(prefix) #write all the output to this file
    
    tntstring += "\ntaxname =;" #les noms de taxons ne sont pas remplacés par des chiffres
    
    if analysis == "bandb":
        tntstring += "\nienum;" #branch and bound
        
    elif analysis == "heuristic":
        tntstring += "\nhold 10000;" #max trees
        tntstring += "\nmult=replic {};".format(str(nrep)) #heuristic
        tntstring += "\nbbreak=tbr;" #swapping method
        tntstring += "\ncollapse [; collapse 1;" #collapse null branches
        
    tntstring += "\nexport - {}.tre;".format(prefix) #export the MPRs and the consensus tree to a file
    tntstring += "\nquit"

    if count_zeroweights > 0:
        print("Warning: {} zero weight triplets deleted due to rounding".format(count_zeroweights))

    with open(prefix+".tnt", "w") as tnt_file:
        tnt_file.write(tntstring)


def triplet_tmc_file(triplet_dict, character_dict, taxa_dict, taxa_convers, prefix=None, weighting="FW"): #construction d'un fichier pour TMC à partir d'une liste de triplets + FW
    """
        Produit un fichier lisible par triplet maxcut (TMC) à partir d'un
        fichier contenant des arbres en format newick plus un fichier contenant
        la correspondance des nombres
    """
    if prefix==None:
        root = Tk()
        prefix =  filedialog.asksaveasfilename(title = "Save three-taxon analysis TMC file", 
                                                filetypes = (("text files","*.txt"),("all files","*.*")))
        root.withdraw()

    #fichier translation taxa names to numbers
    with open(prefix+"_taxa.txt", "w") as TMC_file:
        for convers_line in taxa_convers:
            TMC_file.write(convers_line+"\n") #lignes de conversion taxons/nombres
        TMC_file.write("\n")
    
    #string file
    tmcstring = ''
    for triplet1 in triplet_dict.keys():
        
        in_tax_list = list(triplet1.in_taxa)
        (out_taxa_str,) = triplet1.out_taxa

        #if weights are fractions
        if weighting == "FW" or weighting == "FWNL" or weighting == "MW":
            if character_dict:
                FWdec = len(str(len(character_dict)))
            else: 
                FWdec = 2

            tmcstring += taxa_dict[in_tax_list[0]]+","+taxa_dict[in_tax_list[1]]+"|"+taxa_dict[out_taxa_str]+"[:"+str(round((float(triplet1.FW)),8-FWdec))+"] "
        
        #if weights are integers
        else:
            tmcstring += taxa_dict[in_tax_list[0]]+","+taxa_dict[in_tax_list[1]]+"|"+taxa_dict[out_taxa_str]+"[:"+str(int(triplet1.FW))+"] "

    #core TMC file
    with open(prefix+".tmc", "w") as TMC_file:
        TMC_file.write(tmcstring)
        

def triplet_wqfm_file(triplet_dict, prefix, multiplier=1000000, weighting="FW"): #construction d'un fichier pour TMC à partir d'une liste de triplets + FW
    """
        Produit un fichier lisible par wqfm à partir d'une liste de triplets pondérés
    """

    #core wqfm file
    with open(prefix+".wqfm", "w") as wqfm_file:
        for triplet1, FW in triplet_dict.items():
                        
            #if weights are fractions
            if weighting == "FW" or weighting == "FWNL" or weighting == "MW":
                wqfm_file.write("(("+list(triplet1.out_taxa)[0]+",root),("+list(triplet1.in_taxa)[0]+","+list(triplet1.in_taxa)[1]+")); "+str(round(float(FW)*int(multiplier)))+"\n")
            
            #if weights are integers
            else:
                wqfm_file.write("(("+list(triplet1.out_taxa)[0]+",root),("+list(triplet1.in_taxa)[0]+","+list(triplet1.in_taxa)[1]+")); "+str(round(FW))+"\n")


def triplet_to_file(triplet_dict, character_dict, prefix, analysis="heuristic", 
                    nreplicates=1000, logfileb=True, software="paup", 
                    weighting="FW", wmul=1000):
    """
    Prend un dico de triplets et écrit un fichier tnt/paup...

    Parameters
    ----------
    triplet_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    print("Computing {} file.".format(software))
    
    start = time.time()
    
    if analysis == "auto":
        if character_dict:
            if len(taxa_list_extraction(character_dict)) > 15:
                analysis = "heuristic"
                print("Analysis automatically set to heuristic")
            else:
                analysis = "bandb"
                print("Analysis automatically set to branch and bound")
        else:
            analysis = "bandb"
            print("Analysis automatically set to branch and bound")

    #écriture du fichier en fonction du format demandé
    if software == "paup":
        triplet_nexus_file(triplet_dict, character_dict, weighting, analysis, prefix, 
                           nrep=nreplicates, logfile=logfileb)
    elif software == "tmc":
        taxa_list = taxa_list_extraction(character_dict)
        taxa_dict, taxa_convers = taxa_to_numbers(taxa_list)
        triplet_tmc_file(triplet_dict, character_dict, taxa_dict, taxa_convers, 
                         prefix=None, weighting="FW")
    elif software == "tnt":
        triplet_tnt_file(triplet_dict, character_dict, weighting, analysis, prefix=prefix, 
                         nrep=nreplicates,logfile=logfileb)
    elif software == "wqfm":
        triplet_wqfm_file(triplet_dict, prefix=prefix, multiplier=wmul, weighting=weighting)
    
    end = time.time()
    time_cptr = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print("elapsed time (output file computation): {}".format(time_cptr))


def convert(infile, infiletype, prefix, parallel="auto", 
            weighting="FW", analysis="heuristic", taxa_replacement=False, 
            taxa_replacement_file=False, nreplicates=1000, logfileb=True, 
            software="paup", verbose=True, wmul=1000):
    """
    Prend soit un chemin vers un fichier d'arbres, soit un character_dict,
    soit un triplet dict
    Créé une matrice de 3is au format désiré

    Parameters
    ----------
    infile : TYPE
        DESCRIPTION.
    infiletype : TYPE
        DESCRIPTION.
    prefix : TYPE
        DESCRIPTION.
    parallel : TYPE, optional
        DESCRIPTION. The default is "auto".
    weighting : TYPE, optional
        DESCRIPTION. The default is "FW".
    analysis : TYPE, optional
        DESCRIPTION. The default is "heuristic".
    taxa_replacement : TYPE, optional
        DESCRIPTION. The default is False.
    taxa_replacement_file : TYPE, optional
        DESCRIPTION. The default is False.
    nreplicates : TYPE, optional
        DESCRIPTION. The default is 1000.
    logfileb : TYPE, optional
        DESCRIPTION. The default is True.
    software : TYPE, optional
        DESCRIPTION. The default is "paup".
    verbose : TYPE, optional
        DESCRIPTION. The default is True.
    wmul : TYPE, optional
        DESCRIPTION. The default is 1000.

    Returns
    -------
    None.

    """

    if infiletype == "triplets":
        triplet_dict = triplet_extraction(infile, taxa_replacement_file=taxa_replacement)   
        character_dict = False         
        
    elif infiletype == "trees":
        if type(infile) == str:
            character_dict = character_extraction(infile, taxa_replacement=False)
        
        else:
            character_dict = infile
            
        triplet_dict = main_tripdec(infile, prefix, taxa_replacement, 
                                      weighting, parallel, verbose)

    triplet_to_file(triplet_dict, character_dict, prefix, analysis=analysis, 
                    nreplicates=nreplicates, logfileb=True, software=software, 
                    weighting=weighting, wmul=1000)

        
def nexus_mulcol(infile, prefix, wmultiplier=100):
    """
    Prend un fichier nexus de triplets avec des poids (format PAUP: wts)
    Renvoie un fichier nexus sans poids mais avec les colonnes multipliées par les poids
    La variable wmultiplier correspond à la valeur à laquelle on doit multiplier les poids
    avant d'opérer l'écriture du nexus (les poids sont le plus souvent des float > 1)

    """

    #extraction de la matrice et des poids
    with open(infile, "r") as nexw:
        matrixcptr = False
        wtscptr = False
        matrix = []
        matrixnames = []
        for line in nexw:
    
            #détection de la matrice
            if matrixcptr == False: 
                if "Matrix" in line:
                    matrixcptr = True
                    
            else: 
                if ";" in line:
                    matrixcptr = False
                else:
                    linesplit = line.split(" ")
                    if len(linesplit) == 2: #si la ligne n'est pas vide
                        matrixnames.append(linesplit[0]) #ajout des lignes dans matrix
                        matrix.append(linesplit[1]) #ajout des lignes dans matrix
    
            #détection de la ligne des poids
            if wtscptr == False:
                if "wts" in line:
                    wtscptr = True
            
            else:
                if ";" in line:
                    wtscptr = False
                else:
                    wtslist = [round(float(e.split(":")[0])*wmultiplier) for e in line.split(", ")]

    #s'il y a des poids à zéro
    zerowts = wtslist.count(0) #compteur de poids à zéro
    if zerowts > 0:
        print("WARNING: {} too low weights raised to 1 du to PAUP max weight value limitation".format(zerowts))
        nwtslist = [i if i > 0 else 1 for i in wtslist]
        wtslist = nwtslist

    #écriture du nouveau fichier
    with open(prefix+".nex", "w") as nexus_file:
        nexus_file.write("#NEXUS")
        nexus_file.write("\nbegin data;")
        nexus_file.write("\nDimensions ntax={} nchar={};".format(len(matrixnames),sum(wtslist)))
        nexus_file.write("\nFormat symbols=\"0 1\" missing=?;")
        nexus_file.write("\n")
        nexus_file.write("\nMatrix")
        nexus_file.write("\n")
        
        for i in range(len(matrixnames)): #pour chaque ligne (i correspond au numéro de ligne)
            nexus_file.write("\n"+matrixnames[i]+" ") #début de ligne avec taxon et espace
            newline = "".join([char*count for char,count in zip(matrix[i],wtslist)]) #multiplication des colonnes ligne par ligne
            nexus_file.write(newline) #écriture de la ligne
            
        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")
        nexus_file.write("\n")
        nexus_file.write("\nBegin Paup;")
        nexus_file.write("\noutgroup root /only;")
        nexus_file.write("\nhsearch addseq=random nreps=1000;")
        nexus_file.write("\nroottrees;")
        nexus_file.write("\nsavetrees /file={}.tre format=newick".format(prefix))
        nexus_file.write("\n;")
        nexus_file.write("\n")
        nexus_file.write("\nend;")


# def tree_lisbeth_file(fichier, character_dict): 
#     """
#         permet de générer un fichier 3ia lisible par lisbeth à partir d'un 
#         character dictionary
#     """
#     taxa_dict, taxa_txt_bloc = taxa_to_lisbeth_numbers(taxa_list_extraction(character_dict))
#     for character in character_dict.keys():     #remplacement des noms de taxons par des identifiants
#         for leaf in character.iter_leaves():
#             leaf.name = taxa_dict[leaf.name]

#     #génération du fichier 3ia (Lisbeth)
#     with open(fichier, "w") as tta_file:
#         tta_file.write(";")
#         tta_file.write("\nTaxa")
#         for line in taxa_txt_bloc:
#             tta_file.write("\n" + line)
        
#         tta_file.write("\n;")
#         tta_file.write("\nCharacters")
        
#         for character, char_nb in character_dict.items():
#             char_line = character.write(format=9)
#             char_line = char_line.replace(";", "")
#             char_line = char_line.replace(",", " ")
#             tta_file.write("\n[{}] {}".format(char_nb, char_line))

#         tta_file.write("\n;")

#     for character in character_dict.keys(): #on remet le character_dict à la normale (les noms des taxons des arbres)
#         for leaf in character.iter_leaves():
#             for taxa_name1, symbol in taxa_dict.items():
#                 if leaf.name == symbol:
#                     leaf.name = taxa_name1

#     return taxa_dict, taxa_txt_bloc, character_dict


def agatta_analysis(file_path, software_path, software="tnt", 
                    taxa_replacement=False, method="No", weighting="FW", 
                    parallel="auto", prefix="agatta_out", analysis="bandb", 
                    nrep=1000, rosette=False, chartest=False, ri=False, 
                    consensus=False, pdf_file=False, verbose=True):
    
    print("Starting analysis.")
    
    with open(prefix+".log", "w") as log_file:
        log_file.write("Analysis parameter\n\n")
        log_file.write("Current date and time : ")
        now = datetime.datetime.now()
        log_file.write(now.strftime("%Y-%m-%d %H:%M:%S"))   
        log_file.write("\n")
        log_file.write("Input file path: "+file_path+"\n")
        log_file.write("Prefix: "+prefix+"\n")
        log_file.write("Taxa replacement:"+str(taxa_replacement)+"\n")
        log_file.write("Method: "+method+"\n")
        
        if rosette:
            log_file.write("Standardisation: yes\n")
            log_file.write("Table path for standardisation: "+rosette+"\n")
        else:
            log_file.write("Standardisation: no\n")
            
        log_file.write("Software used: "+software+"\n")   
        
        if software_path:
            log_file.write("Software path: "+software_path+"\n")
        else:
            log_file.write("Software path: no file path\n")
            
        log_file.write("Weighting triplets: "+weighting+"\n")
        log_file.write("Parallelisation: "+parallel+"\n")            
            
        if analysis == "heuristic":
            log_file.write("Analysis: heuristic with {} replicates\n".format(str(nrep)))
        else:
            log_file.write("Analysis: branch and bound\n")
            
        log_file.write("Character states test: "+str(chartest)+"\n")
        
        if pdf_file:
            log_file.write("Character states test pdf results location: "+str(pdf_file)+"\n")
        else:
            log_file.write("Character states test: no pdf required by user\n")
            
        log_file.write("Retention index: "+str(ri)+"\n")
        log_file.write("Consensus: "+str(consensus)+"\n")
        log_file.write("Verbose: "+str(verbose)+"\n")

    #reconnaitre si input est hmatrix ou treelist (.hmatrix)
     #si hmatrix, convertir en treelist
    if "hmatrix" in file_path:
        
        f_path = os.path.split(file_path)[0]
        
        character_dict, poly, error_d = hmatrix(file_path,
                                               f_path+"agatta_character_trees",
                                               prefix, verbose)
    
        file_path = os.path.split(file_path)[0]+"agatta_character_trees.tre"
        
    else:
        character_dict = character_extraction(file_path, taxa_replacement)

    #détecter si standardisation nécessaire (option + rosette)
    if rosette:
        
        #save character_dict
        character_dict = phylo_to_areagram(character_dict, 
                          rosette, 
                          prefix,
                          verbose=verbose)
     
    #détecter si répétitions et les supprimer (le signaler à l'utilisateur)
    if rep_detector(character_dict):
        character_dict = del_replications_forest(character_dict,
                         method=method,
                         prefix=prefix,
                         verbose=verbose)

    prefix_path = os.getcwd() + "/" + prefix

    if software == "agatta":
        triplet_dict = main_tripdec(character_dict, prefix, taxa_replacement, 
                                    weighting, parallel, verbose)

        optimal_score, results_dict = bandb(list(taxa_list_extraction(character_dict)), triplet_dict)
        cladogram_dict = {t : 1 for t in results_dict}
        
    #décomposer treelist en triplets et sauvegarde du fichier nex/tnt + log triplets
    else:
        convert(character_dict, 
                "characters", 
                prefix, 
                parallel, 
                weighting, 
                analysis, 
                taxa_replacement, 
                False, 
                nrep, 
                True, #logfile
                software, 
                verbose)
                
        #analyser le fichier via paup, tnt, ou agatta
        if software == "paup":
            prefix_end = prefix_path+".nex"
            search_pipeline(prefix_end, 
                            software_path,#exception if software_path == False
                            "paup")
        
        elif software == "tnt":
            #prefix_end = prefix_path+".tnt"
            search_pipeline(prefix+".tnt", 
                            software_path, 
                            "tnt")
    

        #récupérer le fichier résultats avec les arbres et supprimer profil
        results_dict = dict()
        i = 1
        with open(prefix+".tre","r") as result_file:
            for ln in result_file:
                if ln.startswith("("):
                    tree = Tree(ln.replace(" ", ""))
                    results_dict[tree] = i
                    i += 1
        
        cladogram_dict = dict()

        for cladogram, i in results_dict.items():
            
            root = cladogram.search_nodes(name="root")[0]
    
            root.detach() # = C.remove_child(J)
            
            cladogram_dict[cladogram] = i
    
    ################################################OPTIONS
    
    #générer un consensus
    if consensus:
        
        #strict consensus
        if consensus == "strict":
        
            constrict(list(cladogram_dict.keys()), prefix)
        
        #reduced cladistic consensus
        elif consensus == "rcc":
        
            rcc(list(cladogram_dict.keys()), prefix)
        
    #test des états de caractères sur le consensus strict
    if chartest:
        
        # folder = os.path.join(f_path, "character_states_test")        
                
        # try:
        #     os.mkdir(folder)
        # except OSError:
        #     print ("Creation of the directory %s failed" % folder)
        # else:
        #     print ("Successfully created the directory %s" % folder)
        
        character_states_test(character_dict, 
                              {constrict(list(cladogram_dict.keys())):1},
                              prefix,
                              pdf_file)
        
    #calcul des indices de rétention
    if ri:
        RI(list(cladogram_dict.keys())[0],
           character_dict,
           weighting,
           prefix)

    print('The analysis ended successfully')
