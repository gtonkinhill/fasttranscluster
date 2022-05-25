from fasttranscluster.pairsnp import run_pairsnp
import tempfile
import os

def test_pairsnp(datafolder):
    # test the interface with the c++ version of pairsnp
    fd, path = tempfile.mkstemp()
    try:
        os.close(fd)
        sample_to_index = {'seq1':0, 'seq2':1, 'seq3':2, 'seq4':3, 'seq4':5, 'seq5':6}
        distances = run_pairsnp(datafolder + 'ambig.aln', 10, path, sample_to_index, ncpu=1)
        print(distances)
        distances = [d[2] for d in distances]
    finally:
        os.remove(path)

    assert distances == [0, 2, 1, 1, 2, 2, 2, 3, 3, 0]

    return