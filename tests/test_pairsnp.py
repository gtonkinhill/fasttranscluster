from fasttranscluster.pairsnp import run_pairsnp
import tempfile
import os

def test_pairsnp(datafolder):
    # test the interface with the c++ version of pairsnp
    fd, path = tempfile.mkstemp()
    try:
        os.close(fd)
        distances = run_pairsnp(datafolder + 'ambig.aln', 10, path, ncpu=1)
    finally:
        os.remove(path)

    assert distances == [[0, 1, 0], [0, 2, 2], [0, 3, 1], [0, 4, 1], [1, 2, 2], [1, 3, 2], [1, 4, 2], [2, 3, 3], [2, 4, 3], [3, 4, 0]]

    return