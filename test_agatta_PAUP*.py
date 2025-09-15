# -*- coding: utf-8 -*-
"""

    Agatta: Three-item analysis Python package
    Contact: Valentin Rineau - valentin.rineau@sorbonne-universite.fr

    Agatta is a set of tools in the cladistic framework to perform
    three-item analysis and associated operations in cladistics in python.

    https://github.com/vrineau/agatta

    This code is under license GNU GPLv3

"""

import os
import shutil
from ete3 import Tree
from tqdm import tqdm
from fractions import Fraction
from pathlib import Path
from agatta.main import coremain

#The following tests are using a binary version of PAUP*

def test_gill_morpho(tmp_path):
    
    # temp directory 
    cwd_backup = os.getcwd()
    os.chdir(tmp_path)
    
    # test file
    data_dir = Path(__file__).parent / "data"
    source_file = data_dir / "gill_morpho.hmatrix"

    # Copy to tmp_path
    test_file = tmp_path / "gill_morpho.hmatrix"
    test_file.write_text(source_file.read_text())

    # PAUP path, work also on Github Actions
    paup_path = shutil.which("paup")  
    if paup_path is None:
        paup_path = "/home/valentin/paup/paup4a168_ubuntu64"
        
    # prefix
    prefix = "gill_morpho_results"

    # AGATTA arguments
    arguments = {'--analysis': 'heuristic',
                 '--chardec': False,
                 '--chartest': False,
                 '--consensus': 'strict',
                 '--detailed_tripdec': False,
                 '--directory': './',
                 '--filetype': 'trees',
                 '--help': False,
                 '--index': 'ri',
                 '--log': False,
                 '--nsupport': False,
                 '--parallel': 'auto',
                 '--pdf': False,
                 '--prefix': prefix,
                 '--repetitions': 'TMS',
                 '--replicates': '1000',
                 '--ri': False,
                 '--rnri_codes': None,
                 '--rnri_rescaling': False,
                 '--rnri_totaltree': False,
                 '--rosetta': None,
                 '--showtaxanames': False,
                 '--softpath': paup_path,
                 '--software': 'paup',
                 '--taxarep1': None,
                 '--taxarep2': None,
                 '--version': False,
                 '--weighting': 'FW',
                 '-s': False,
                 '-v': False,
                 '<command>': None,
                 '<file>': [str(test_file)],
                 'analysis': True,
                 'chartest': False,
                 'consensus': False,
                 'convert': False,
                 'describetree': False,
                 'fp': False,
                 'help': False,
                 'hmatrix': False,
                 'nri': False,
                 'standardisation': False,
                 'support': False,
                 'tripdec': False}

    # launch analysis
    try:
        coremain(arguments)
    finally:
        os.chdir(cwd_backup)

    # check output files
    output_file = tmp_path / f"{prefix}.tre"
    assert output_file.exists(), "The output file was not created."

    # # Check content
    # content = output_file.read_text()
    # assert "1 optimal tree found" in content
    # assert "Strict consensus computed" in content
