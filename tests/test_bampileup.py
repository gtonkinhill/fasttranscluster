from fasttranscluster.bampileup import preprocess_pileup
import tempfile
import os
import gzip
import numpy as np

def test_preprocess_pileup(datafolder):
    # A very simple test of the pileup function. 
    # Does not check estimated probabilities are correct.
    fd, path = tempfile.mkstemp()
    try:
        os.close(fd)
        preprocess_pileup([datafolder + 'test.bam'],
                datafolder + 'test_ref.fasta',
                path,
                expected_freq_threshold=None,
                keep_all=True,
                require_both_strands=False,
                min_contig_length=1000,
                max_depth=8000,
                min_mapq=60,
                min_baseq=13,
                quiet=False)
        with gzip.open(path, 'rt') as infile:
            pos = infile.read().strip()
            pos = pos.strip().split(',')
            sample_name = pos[0]
            a = np.array(pos[1:]).astype(float)
            a = np.reshape(a, (int(a.size/4),4), order='C')
    finally:
        os.remove(path)

    assert sample_name == 'test'
    assert np.sum((a>1) | (a<0)) == 0

    return