from fasttranscluster.__main__ import main
from fasttranscluster.bampileup import preprocess_pileup
import tempfile
import os, sys
import numpy as np

def test_ftclust_msa(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        
        # run ftclust
        sys.argv = [
            "", "--msa", datafolder + "/ambig.aln",
            "--dates", datafolder + 'dates_msa.csv', 
            "-o", tmpoutdir, "-K", "10", "-P", "0", "--snp_threshold", "5", "--save_probs"
        ]
        main()

        with open(tmpoutdir + '/transcluster_probabilities.csv', 'r') as infile:
            lines = infile.read()

        # read probability file
        assert "seq1,seq2,0" in lines
        assert "seq3,seq4,3" in lines

    return