from fasttranscluster.__main__ import main
from fasttranscluster.bampileup import preprocess_pileup
import tempfile
import os, sys
import numpy as np

def test_ftclust_pileup(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:

        for b in ['test', 'testB']:
            preprocess_pileup([datafolder + b +'.bam'],
                    datafolder + 'test_ref.fasta',
                    tmpoutdir + '/' + b + '.csv.gz',
                    expected_freq_threshold=0.1,
                    keep_all=False,
                    require_both_strands=False,
                    min_contig_length=1000,
                    max_depth=8000,
                    min_mapq=60,
                    min_baseq=13,
                    quiet=False)


        # run ftclust
        sys.argv = [
            "", "--pileup", tmpoutdir + "/test.csv.gz", tmpoutdir + "/testB.csv.gz",
            "--dates", datafolder + 'test_dates.csv', 
            "-o", tmpoutdir, "-K", "10", "-P", "0", "--snp_threshold", "99999", "--save_probs"
        ]
        main()

        with open(tmpoutdir + '/transcluster_probabilities.csv', 'r') as infile:
            lines = infile.read()

        # read probability file
        assert "test,testB,2" in lines

    return