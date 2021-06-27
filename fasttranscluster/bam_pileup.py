import os
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import pysam
import pysamstats
import pyfastx
import numpy as np
from fasttranscluster.dirichlet_multinomial import find_dirichlet_priors


def pileup_dist(alnfiles,
                reference,
                outputfile,
                expected_freq_threshold=0.01,
                min_hit_cov=3,
                require_both_strands=True,
                min_contig_length=1000,
                max_depth=8000,
                quiet=False):

    # load pileups as output by pysamstats
    if not quiet:
        print('Loading bams...')

    # load reference contig lengths
    contig_lengths = []
    total_length = 0
    for name, seq in pyfastx.Fasta(reference, build_index=False):
        l = len(seq)
        if l < min_contig_length: continue
        contig_lengths.append((name, l))
        total_length += l

    attr = [(0, 'A_fwd', 'A_rev'), (1, 'C_fwd', 'C_rev'), (2, 'G_fwd', 'G_rev'), (3, 'T_fwd', 'T_rev')]

    # 0 (no evidence) /  1 (found) / 9 (unknown)
    with open(outputfile, 'ab') as outfile:
        for alnfile in alnfiles:
            prefix,ext = os.path.splitext(os.path.basename(alnfile))
            if ext=='sam':
                samfile = pysam.AlignmentFile(alnfile, "r")
            elif ext=='sam':
                samfile = pysam.AlignmentFile(alnfile, "rb")
            elif ext=='cram':
                samfile = pysam.AlignmentFile(alnfile, "rc")
            else:
                raise ValueError('File extension is not supported!')

            all_counts = np.zeros((total_length, 4), dtype=int)

            current = 0
            for name, l in contig_lengths:
                print(name, l)
                pile = pysamstats.load_variation_strand(samfile, 
                    chrom=name, 
                    start=1,
                    end=l,
                    truncate=False,
                    fafile=reference,
                    max_depth=max_depth)
                for i, f, r in attr:
                    counts = pile[f] + pile[r]
                    if require_both_strands:
                        counts[pile[f]<=0] = 0
                        counts[pile[r]<=0] = 0
                    all_counts[current+pile['pos'], i] = counts
                current += l

            alphas = find_dirichlet_priors(all_counts)
            a0 = np.sum(alphas)
            denom = a0 + np.sum(all_counts, 1)

            all_counts = (all_counts+a0)/denom

            # Filter out those below the threshold
            all_counts[all_counts<threshold] = 0
            
            # save output to file
            outfile.write(prefix + ',')
            np.savetxt(outfile, all_counts, sep=',', newline=',',format='%0.3f')
            outfile.write('\n')

    # # calculate pairwise distances
    # if not quiet:
    #     print('calculating distances...')
    # distances = []
    # for i in tqdm(range(len(pileups))):
    #     for j in range(i + 1, len(pileups)):
    #         d = 0
    #         for chrom in variants[i]:
    #             d += len(index) * MAX_LEN - np.unpackbits(
    #                 (variants[i][0] & variants[j][0])
    #                 | (variants[i][1] & variants[j][1])
    #                 | (variants[i][2] & variants[j][2])
    #                 | (variants[i][3] & variants[j][3])).sum()
    #         if d < mid_dist:
    #             distances.append((i, j, int(d)))

    return
