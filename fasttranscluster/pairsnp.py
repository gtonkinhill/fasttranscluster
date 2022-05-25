import os
import subprocess


def run_pairsnp(msa, snp_threshold, outputfile, sample_to_index, ncpu=1):
    # runs pairsnp and reads result into a numpy array

    cmd = "pairsnp"
    cmd += " -d " + str(snp_threshold)
    cmd += " -t " + str(ncpu)
    cmd += " -s"
    cmd += " " + msa
    cmd += " > " + outputfile

    subprocess.run(cmd, shell=True, check=True)

    # load output
    distances = []
    with open(outputfile, "r") as infile:
        for line in infile:
            line = line.strip().split()
            distances.append(
                [sample_to_index[line[0]], sample_to_index[line[1]], int(line[2])]
            )

    return distances
