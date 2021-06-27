import os
import datatable as dt
import numpy as np
from collections import defaultdict
from tqdm import tqdm

MAX_LEN = 1000000


def pileup_dist(pileups,
                ignoreN=True,
                min_ref_support=1,
                mid_dist=9999999999999,
                is_min=True,
                quiet=False):

    # load pileups as output by pysamstats
    if not quiet:
        print('Loading pileups...')

    sample_names = []
    variants = []
    ignore_pos = []
    index = {}
    index_count = 0

    for i, f in tqdm(enumerate(pileups)):
        sample_names.append(
            os.path.splitext(os.path.basename(f).replace('.gz', ''))[0])

        variants.append([])
        # variants.append(defaultdict(list))
        d = dt.fread(f)
        d[:,
          dt.update(As=(dt.f.A_fwd >= 1) & (dt.f.A_rev >= 1)
                    & ((dt.f.A_fwd + dt.f.A_rev) >= 1))]
        d[:,
          dt.update(Cs=(dt.f.C_fwd >= 1) & (dt.f.C_rev >= 1)
                    & ((dt.f.C_fwd + dt.f.C_rev) >= 1))]
        d[:,
          dt.update(Gs=(dt.f.G_fwd >= 1) & (dt.f.G_rev >= 1)
                    & ((dt.f.G_fwd + dt.f.G_rev) >= 1))]
        d[:,
          dt.update(Ts=(dt.f.T_fwd >= 1) & (dt.f.T_rev >= 1)
                    & ((dt.f.T_fwd + dt.f.T_rev) >= 1))]
        d[:, dt.update(has_entry=(dt.f.As | dt.f.Cs | dt.f.Gs | dt.f.Ts))]

        chroms = dt.unique(d['chrom']).to_numpy().flatten()
        d = d[dt.f.has_entry == 1, :]

        for c in chroms:
            if c not in index:
                index[c] = index_count
                index_count += 1

        ind = np.zeros(d.shape[0], dtype=int)
        for j, chrom in enumerate(d['chrom'].to_list()[0]):
            ind[j] = index[chrom]

        if i == 0:
            trim_length = MAX_LEN * len(index)

        array_index = ind * MAX_LEN + d['pos'].to_numpy().flatten()
        temp_array = np.ones(MAX_LEN * len(index), dtype=bool)
        temp_array[array_index] = d['As'].to_numpy().flatten()
        variants[i].append(np.packbits(temp_array[:trim_length]))
        temp_array[array_index] = 1

        temp_array[array_index] = d['Cs'].to_numpy().flatten()
        variants[i].append(np.packbits(temp_array[:trim_length]))
        temp_array[array_index] = 1

        temp_array[array_index] = d['Gs'].to_numpy().flatten()
        variants[i].append(np.packbits(temp_array[:trim_length]))
        temp_array[array_index] = 1

        temp_array[array_index] = d['Ts'].to_numpy().flatten()
        variants[i].append(np.packbits(temp_array[:trim_length]))
        temp_array[array_index] = 1

    # calculate pairwise distances
    if not quiet:
        print('calculating distances...')
    distances = []
    for i in tqdm(range(len(pileups))):
        for j in range(i + 1, len(pileups)):
            d = 0
            for chrom in variants[i]:
                d += trim_length - np.unpackbits(
                    (variants[i][0] & variants[j][0])
                    | (variants[i][1] & variants[j][1])
                    | (variants[i][2] & variants[j][2])
                    | (variants[i][3] & variants[j][3])).sum()
            if d < mid_dist:
                distances.append((i, j, int(d)))

    return (distances)
