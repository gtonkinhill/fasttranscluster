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
        print(i)
        sample_names.append(
            os.path.splitext(os.path.basename(f).replace('.gz', ''))[0])
        variants.append(defaultdict(lambda: set('R')))
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
        poss = d['pos'].to_numpy().flatten()
        chroms = d['chrom'].to_numpy().flatten()
        ignore_pos.append(
            set(zip(chroms[has_entry == False], poss[has_entry == False])))

        ref = d['ref'].to_numpy().flatten()
        npos = len(ref)
        for j in np.arange(npos)[(ref != 'A') & As]:
            variants[i][(chroms[j], poss[j])].add('A')
        for j in np.arange(npos)[(ref != 'C') & Cs]:
            variants[i][(chroms[j], poss[j])].add('C')
        for j in np.arange(npos)[(ref != 'G') & Gs]:
            variants[i][(chroms[j], poss[j])].add('G')
        for j in np.arange(npos)[(ref != 'T') & Ts]:
            variants[i][(chroms[j], poss[j])].add('T')

        # remove default reference from positions where its missing
        for j in np.arange(npos)[(ref == 'A') & (As == False)]:
            variants[i][(chroms[j], poss[j])].remove('R')
        for j in np.arange(npos)[(ref == 'C') & (Cs == False)]:
            variants[i][(chroms[j], poss[j])].remove('R')
        for j in np.arange(npos)[(ref == 'G') & (Gs == False)]:
            variants[i][(chroms[j], poss[j])].remove('R')
        for j in np.arange(npos)[(ref == 'T') & (Ts == False)]:
            variants[i][(chroms[j], poss[j])].remove('R')

    # calculate pairwise distances
    if not quiet:
        print('calculating distances...')
    distances = []
    for i in range(len(pileups)):
        print(i)
        for j in range(i + 1, len(pileups)):
            d = 0
            for pos in (variants[i].keys()
                        | variants[j].keys()) - ignore_pos[i] - ignore_pos[j]:
                if len(variants[i][pos].intersection(variants[j][pos])) <= 0:
                    d += 1
            if d < mid_dist:
                distances.append((i, j, d))
                distances.append((j, i, d))

    return (distances)
