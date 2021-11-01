# AGATTA - Three-item analysis python package

Three-item analysis python package. 3is matrix computation.

## Installation:
	This package can be installed using pip install filename.wh


## Usage:    
    agatta tripdec <file> [-s -i --prefix=<file>]
    agatta mulcol <file> [-s --prefix=<file> --multiplier=<int>]
    agatta ri <file> <file> [-s -i -j --prefix=<file>]
    agatta itri <file> <file> [-s -i -j --prefix=<file>]
    agatta itrisym <file> <file> [-s -i -j --itrimethod=<type> --prefix=<file>]
    agatta chartest  <file> <file> [-s -i -j --directory=<dir>]
    agatta convert (--paup | --tmc | --tnt) <file> [-s -i -v -p=<type> --log --weighting=<type>  --nrep=<int> --analysis=<str> --prefix=<file>]
    agatta fp <file> [-s -v -i --method=<type> --prefix=<file>]
    agatta constrict <file> [-s -i --prefix=<file>]
    agatta rcc <file> [-s -i --prefix=<file>]
    agatta describetree <file> [-s --prefix=<file>]
    agatta areaconvert <file> <file> [-s -v]
    agatta hmatrix <file> [-s -v --chardec --method=<type> --prefix=<file>]
    agatta -h | --help

## Options:
    -h --help            show this
    -v                   verbose mode
    -s                   silent mode
    -i                   replace taxa names by id from taxabloc for the first file
    -j                   replace taxa names by id from taxabloc for the second file
    -p=<type>            fast tree decomposition using parallelisation (choose auto for automatic detection of the number of cores, or choose a specific number of cores)
    --prefix=<file>      prefix of the output file [default: agatta_out]
    --directory=<dir>    directory for output chartest files [default: ./]
    --paup               compute a nexus file for PAUP from a list of character trees 
    --tmc                compute a TMC file from a list of character trees with a specified maximal weight
    --tnt                compute a nexus file for TNT from a list of character trees 
    --weighting=<type>   specify the type of triplet weighting for the three-item analysis (FW,FWNL,UW,MW,AW,NW) [default: FW]
    --method=<type>      specify the method used to produce subtrees (VR, RZB, No) [default: VR]
    --itrimethod=<type>  specify the method used to compute symmetric ITRI (sum, product) [default: sum]
    --multiplier=<int>   specifies an integer by which the weight of each column in the PAUP file is multiplied [default: 1]
    --log                add a line to save a log file
    --chardec            decompose character trees into components from a cao matrix
    --analysis=<str>     type of analysis, heuristic or bandb [default: heuristic]
    --nrep=<int>         number of replicates for an analysis [default: 1000]

 
## Exhaustive list of available commands:

tripdec		decompose a character set into a list of triplets with their fractional weights")
		
		Input(one file requested): 
		characters file. 

		Output (file name requested): 
		list of triplets with their fractional weights


mulcol        convert a matrix with weights to a matrix with column repeated (number of repetition equals rounded weight * multiplier) instead of weights.

		Input(one file requested): 
		nexus file generated by agatta or lisbeth. 

		Output (file name requested): 
		nexus file

ri   	 	calculate the retention index of a character set relatively to a cladogram")
		
		Input(two files requested): 
		(1) cladogram file
		(2) characters file
		
		Output (file name requested): 
		list of RI by character

itri		calculate the inter-tree retention index between two trees (BEWARE: ITRI is not symmetrical, to
		obtain a (symmetric) distance measure, take the mean of ITRI 1-2 and ITRI 2-1

		Input(two files requested): 
		(1) tree file with one file
		(2) tree file with one file (reference tree)

		Output:
		values of resolving power (%), artefactual resolution (%) and efficiency (%)

itrisym   triplet distance, i.e. symmetric version of the ITRI (mean of ITRI 1-2 and ITRI 2-1)

		Input(two files requested): 
		(1) tree file with one file
		(2) tree file with one file

		Output:
		values of resolving power (%), artefactual resolution (%) and efficiency (%)
		
chartest   	for a set of character states, compute a pdf file showing reject or acceptation of each character state on a cladogram
		
		Input(two files requested): 
		(1) cladogram file, 
		(2) characters file. 

		Output (folder requested): 
		(3) one pdf file for each character tree in (2)
		(4) list of acceptation/reject by character state and 
		position of the state on the tree (clades labels are on the pdfs)

convert --paup  compute a nexus file for PAUP from a list of character trees
		Input(one file requested): characters file. Output (prefix requested): PAUP nexus file


convert --tmc   compute a Triplet Maxcut file from a list of character trees
		Input(one file requested): characters file. Output (prefix requested): TMC file

convert --tnt
		Input(one file requested): characters file. Output (prefix requested): tnt file


fp     compute a free-paralogy subtree analysis
		Input(one file requested): list of newick trees. Output: list of newick trees without repetitions

constrict          compute a strict consensus from a set of trees
		Input(one file requested): tree file. Output: tree file

rcc                compute a reduced cladistic consensus from a set of trees
		Input(one file requested): tree file. Output: tree file

describetree show various informations on trees
		Input(one file requested): tree file. Output: data file
 
 areaconvert convert a phylogram in areagram using a conversion table (list of taxa:area)
 		Input(two files requested): tree file and conversion table file. Output: tree file


## Examples of use:

    agatta tripdec test2.3ia

    agatta convert --paup --weighting=FW test2.3ia -i 
    agatta convert --paup --weighting=FWNL test2.3ia    
    agatta convert --paup --weighting=UW test2.3ia
    agatta convert --paup --weighting=MW test2.3ia        
    agatta convert --paup --weighting=AW test2.3ia    
    agatta convert --paup --weighting=NW test2.3ia    
    agatta convert --tmc --weighting=FW test2.3ia
    agatta convert --tNT --weighting=FW test2.3ia
    agatta mulcol --multiplier=50 test.nex

	agatta ri test1.tre test2.3ia -j
    agatta itri test1.tre test3.tre
    agatta itrisym test1.tre test3.tre --itrimethod=sum
    agatta itrisym test1.tre test3.tre --itrimethod=product
	
    agatta chartest  test1.tre test2.3ia -j
    agatta describetree test2.3ia

    agatta areaconvert test3.tre areaconvertest.txt -v
    agatta fp testfp.tre -v

	agatta hmatrix flakes_hmatrix_nodoublon.csv -v --prefix=flakes_hmatrix_nodoublon #free-paralogy faite dans la foulée

    agatta --version
    agatta -h

This code is under licence *GNU GPLv3*

















