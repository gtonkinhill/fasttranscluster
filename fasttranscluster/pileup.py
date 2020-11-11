import os
import datatable as dt
import numpy as np
from collections import defaultdict


def pileup_dist(pileups,
                ignoreN=True,
                min_ref_support=3,
                mid_dist=9999999999999,
                is_min=True,
                quiet=False):

    # load pileups as output by pysamstats
    if not quiet:
        print('Loading pileups...')

    sample_names = []
    variants = []
    ignore_pos = []
    for i, f in enumerate(pileups):
        sample_names.append(
            os.path.splitext(os.path.basename(f).replace('.gz', ''))[0])
        variants.append(defaultdict(lambda: np.ones((4, 1000000), dtype=bool)))
        d = dt.fread(f)
        As = ((d['A_fwd'].to_numpy() >= 1) & (d['A_rev'].to_numpy() >= 1) & (
            (d['A_fwd'].to_numpy() + d['A_rev'].to_numpy()) >= 3)).flatten()
        Cs = ((d['C_fwd'].to_numpy() >= 1) & (d['C_rev'].to_numpy() >= 1) & (
            (d['C_fwd'].to_numpy() + d['C_rev'].to_numpy()) >= 3)).flatten()
        Gs = ((d['G_fwd'].to_numpy() >= 1) & (d['G_rev'].to_numpy() >= 1) & (
            (d['G_fwd'].to_numpy() + d['G_rev'].to_numpy()) >= 3)).flatten()
        Ts = ((d['T_fwd'].to_numpy() >= 1) & (d['T_rev'].to_numpy() >= 1) & (
            (d['T_fwd'].to_numpy() + d['T_rev'].to_numpy()) >= 3)).flatten()
        has_entry = (As | Cs | Gs | Ts).flatten()
        poss = d['pos'].to_numpy().flatten() - 1
        chroms = d['chrom'].to_numpy().flatten()
        for chrom in np.transpose(
                d[:, 'chrom', dt.by('chrom')]['chrom'].to_numpy())[0]:
            ind = (chroms == chrom) & has_entry
            variants[i][chrom][0, poss[ind]] = As[ind]
            variants[i][chrom][1, poss[ind]] = Cs[ind]
            variants[i][chrom][2, poss[ind]] = Gs[ind]
            variants[i][chrom][3, poss[ind]] = Ts[ind]

    # calculate pairwise distances
    if not quiet:
        print('calculating distances...')
    distances = []
    for i in range(len(pileups)):
        print(i)
        for j in range(i + 1, len(pileups)):
            d = 0
            for chrom in variants[i]:
                d += np.sum((variants[i][chrom]
                             & variants[j][chrom]).sum(0) == 0)
            if d < mid_dist:
                distances.append((i, j, d))
                distances.append((j, i, d))

    return (distances)
