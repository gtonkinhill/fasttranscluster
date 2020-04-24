import os, sys
import argparse
import pyfastx
from datetime import datetime
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .pairsnp import run_pairsnp
from .transcluster import calculate_trans_prob

from .__init__ import __version__

SECONDS_IN_YEAR = 31556952

def get_options():
    description = 'Runs the pairsnp and transcluster algorithms.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='fasttranscluster')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "--msa",
        dest="msa",
        required=True,
        help="Location of fasta formatted multiple sequence alignment")

    io_opts.add_argument(
        "--dates",
        dest="metadata",
        required=True,
        help="""Location of metadata in csv format. The first column must include the 
        sequence names and the second column must include sampling dates.""")

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    io_opts.add_argument(
        "--tree",
        dest="tree",
        action='store_true',
        default=False,
        help="Toggles the pipeline to build a phylogeny of the initial sequences for each transmission cluster.")

    # transcluster options
    transcluster = parser.add_argument_group('transcluster options')
    transcluster.add_argument(
        "-K",
        "--trans_threshold",
        dest="trans_threshold",
        help=("transmission distance threshold - samples are clustered together" +
            " where the implied number of transmissions k is less than or equal to K " +
            " with a probability of P"),
        type=int,
        default=5)
    
    transcluster.add_argument(
        "-P",
        "--prob_threshold",
        dest="prob_threshold",
        help=("probability threshold - samples are clustered together" +
            " where the implied number of transmissions k is less than or equal to K " +
            " with a probability of P"),
        type=float,
        default=0.8)

    transcluster.add_argument(
        "--clock_rate",
        dest="clock_rate",
        help="clock rate as defined in the transcluster paper (SNPs/genome/year)",
        type=float,
        default=3)

    transcluster.add_argument(
        "--trans_rate",
        dest="trans_rate",
        help="transmission rate as defined in the transcluster paper (transmissions/year)",
        type=float,
        default=52)

    # pairsnp options
    pairsnp = parser.add_argument_group('Pairsnp options')
    pairsnp.add_argument(
        "--snp_threshold",
        dest="snp_threshold",
        help="SNP threshold used to sparsify initial SNP distance matrix",
        type=int,
        default=10)
    
    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    index_to_sample = {}
    sample_to_index = {}
    # get sample names by index from fasta
    for i, seq in enumerate(pyfastx.Fasta(args.msa, build_index=False)):
        index_to_sample[i] = seq[0]
        sample_to_index[seq[0]] = i

    # load metadata
    sample_dates = {}
    with open(args.metadata, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            if line[0] not in sample_to_index: NameError("Missing date in metadata!")
            sample_dates[sample_to_index[line[0]]] = datetime.fromisoformat(line[1]).timestamp()/SECONDS_IN_YEAR

    # run pairsnp
    snp_dist_file = args.output_dir + "pairsnp_sparse_dist.csv"
    sparse_dist = run_pairsnp(msa = args.msa,
        snp_threshold = args.snp_threshold,
        outputfile = snp_dist_file,
        ncpu = args.n_cpu)

    # run transcluster algorithm
    trans_dist = calculate_trans_prob(
        sparse_snp_dist = sparse_dist,
        sample_dates = sample_dates,
        K = args.trans_threshold,
        lamb = args.clock_rate,
        beta = args.trans_rate
    )

    # generate clusters using single linkage algorithm
    sparse_dist_matrix = csr_matrix((data, (row_ind, col_ind)),
                                          shape=(ncentroids, ncentroids))
    sparse_dist_matrix = sparse_dist_matrix <= np.log(args.prob_threshold)
    n_components, labels = connected_components(
        csgraph=sparse_dist_matrix,
        directed=False,
        return_labels=True)

    # write results to file

    # plot results

    return


if __name__ == '__main__':
    main()