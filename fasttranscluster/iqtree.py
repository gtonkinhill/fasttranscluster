import os
import subprocess
import pyfastx


def run_iqtree_index_cases(msa, clusters, dates, outdir, ncpu=1):
    # runs iqtree on the index case for each transmission cluster

    # identify index case
    cluster_indexes = {}
    for i in dates:
        if (clusters[i] not in cluster_indexes) or (
                cluster_indexes[clusters[i]][1] > dates[i][1]):
            cluster_indexes[clusters[i]] = (i, dates[i][1])
    cluster_indexes = set([t[0] for t in cluster_indexes.values()])

    # create fasta file
    index_fasta = outdir + "index_seqs.fasta"
    date_file = outdir + "index_dates.tab"
    with open(index_fasta, 'w') as outfile, open(date_file, 'w') as datefile:
        datefile.write("name\tdate\n")
        for i, seq in enumerate(pyfastx.Fasta(msa, build_index=False)):
            if i in cluster_indexes:
                outfile.write(">cluster_" + str(clusters[i]+1) + "_" + seq[0] +
                              "\n" + seq[1] + "\n")
                datefile.write(seq[0] + "\t" + dates[i][0] + "\n")

    # run iqtree
    cmd = "iqtree2"
    cmd += " -s " + index_fasta
    cmd += " --date " + date_file
    cmd += " -nt " + str(ncpu)
    cmd += " -st DNA"
    cmd += " -quiet"
    cmd += " -redo"

    subprocess.run(cmd, shell=True, check=True)

    return
