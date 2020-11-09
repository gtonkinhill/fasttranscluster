import os
import vcf
import numpy as np
from collections import defaultdict


def vcf_dist(vcfs,
             ignoreN=True,
             min_ref_support=3,
             mid_dist=9999999999999,
             is_min=True,
             quiet=False):

    # load vcfs
    if not quiet:
        print('Loading vcfs...')
    sample_names = []
    variants = []
    ignore_pos = []
    for i, f in enumerate(vcfs):
        sample_names.append(os.path.splitext(os.path.basename(f))[0])
        variants.append(defaultdict(lambda: set('R')))
        ignore_pos.append(set())
        with open(f, 'r') as infile:
            for line in infile:
                if line[0] == '#': continue
                line = line.strip().split()
                dp4 = None
                for k in line[-1].split(';'):
                    temp = k.split('=')
                    if temp[0] == 'DP4':
                        dp4 = [int(c) for c in temp[1].split(',')]
                if dp4 is None:
                    raise ValueError('Missing DP4 attribute in vcf entry!')
                if (dp4[0] > 0) and (dp4[1] > 0) and (
                    (dp4[0] + dp4[1]) >= min_ref_support):
                    variants[i][(line[0], line[1])] |= set(['R', line[4]])
                    # print('ref as well!!')
                else:
                    variants[i][(line[0], line[1])] = set([line[4]])
                    # print('HERE2')

                if ('N' in line[4]) and ignoreN:
                    ignore_pos[i].add((line[0], line[1]))

    variant_positions = [set(v.keys()) for v in variants]

    # calculate pairwise distances
    if not quiet:
        print('calculating distances...')
    distances = []
    for i in range(len(vcfs)):
        print(i)
        for j in range(i + 1, len(vcfs)):
            d = 0
            for pos in (variant_positions[i]
                        | variant_positions[j]) - ignore_pos[i]:
                if len(variants[i][pos].intersection(variants[j][pos])) <= 0:
                    d += 1
            if d < mid_dist:
                distances.append((i, j, d))
                distances.append((j, i, d))

    return (distances)
