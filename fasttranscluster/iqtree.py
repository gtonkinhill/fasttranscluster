import os
import subprocess
import pyfastx


def run_iqtree_index_cases(msa, clusters, dates, outdir, ncpu=1):
    # runs iqtree on the index case for each transmission cluster

    # identify index case
    cluster_indexes = {}
    for i in dates:
        if (clusters[i] not in cluster_indexes) or (
                cluster_indexes[clusters[i]][1] > dates[i][0]):
            cluster_indexes[clusters[i]] = (i, dates[i][0])
    cluster_indexes = set([t[0] for t in cluster_indexes.values()])

    index_fasta = outdir + "index_seqs.fasta"
    with open(index_fasta, 'w') as outfile:
        for i, seq in enumerate(pyfastx.Fasta(msa, build_index=False)):
            if i in cluster_indexes:
                outfile.write(">cluster_" + str(clusters[i]+1) + "_" + seq[0] +
                              "\n" + seq[1] + "\n")

    # run iqtree
    cmd = "iqtree"
    cmd += " -s " + index_fasta
    cmd += " -nt " + str(ncpu)
    cmd += " -st DNA"
    cmd += " -quiet"
    cmd += " -redo"

    subprocess.run(cmd, shell=True, check=True)

    return
