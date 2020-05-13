# fasttranscluster

This is a python implementation of the [transcluster](https://github.com/JamesStimson/transcluster) algorithm. It makes use of [pairsnp](https://github.com/gtonkinhill/pairsnp) as well as a fast summation approach to infer the distribution of transmission probabilities very quickly. It should scale to over 100k genomes.

## Quick Start

`fasttranscluster` requires a multiple sequence alignment as well as a csv file with the dates the samples were obtained. Examples of these can be found in the `example_data` folder in this repository. Currently the default `clock_rate` and `transmission rate` are set to sensible values for sars-cov-2. The clustering can then be controlled by setting the `prob_threshold` and `trans_threshold` parameters.

```
mkdir transcluster_out
fasttranscluster --msa ./example_data/msa.fa --dates ./example_data/meta_data.csv -o ./transcluster_out -K 5
```

The probabilities for each transmission pair that passes the thresholds can be obtained by using the `--save_probs` flag.

## Installation

### Manual

You will first need to install the c++ version of [pairsnp](https://github.com/gtonkinhill/pairsnp-cpp)

Next clone this repsitory and run the setup.py script as

```
git clone https://github.com/gtonkinhill/fasttranscluster
cd fasttranscluster
python setup.py install
```

### Conda

in the future...

## Parameters

```
usage: fasttranscluster [-h] --msa MSA --dates METADATA -o OUTPUT_DIR
                        [--tree {index,mrca}] [-K TRANS_THRESHOLD]
                        [-P PROB_THRESHOLD] [--clock_rate CLOCK_RATE]
                        [--trans_rate TRANS_RATE] [--save_probs]
                        [--snp_threshold SNP_THRESHOLD] [-t N_CPU] [--quiet]
                        [--version]

Runs the pairsnp and transcluster algorithms.

optional arguments:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --quiet               suppress additional output
  --version             show program's version number and exit

Input/output:
  --msa MSA             Location of fasta formatted multiple sequence
                        alignment
  --dates METADATA      Location of metadata in csv format. The first column
                        must include the sequence names and the second column
                        must include sampling dates.
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of an output directory
  --tree {index,mrca}   Toggles the pipeline to build a phylogeny of the
                        initial sequences for each transmission cluster. Can
                        be based on either the first sequnece in each cluster
                        'index' or the MRCA of each cluster 'mrca

transcluster options:
  -K TRANS_THRESHOLD, --trans_threshold TRANS_THRESHOLD
                        transmission distance threshold - samples are
                        clustered together where the implied number of
                        transmissions k is less than or equal to K with a
                        probability of P
  -P PROB_THRESHOLD, --prob_threshold PROB_THRESHOLD
                        probability threshold - samples are clustered together
                        where the implied number of transmissions k is less
                        than or equal to K with a probability of at least P
  --clock_rate CLOCK_RATE
                        clock rate as defined in the transcluster paper
                        (SNPs/genome/year) default=0.9e-3 * 29903
  --trans_rate TRANS_RATE
                        transmission rate as defined in the transcluster paper
                        (transmissions/year) default=73
  --save_probs          write out transmission probabilites (can be a large
                        file)

Pairsnp options:
  --snp_threshold SNP_THRESHOLD
                        SNP threshold used to sparsify initial SNP distance
                        matrix
```

## Citation

Stimson,J., Gardy,J., Mathema,B., Crudu,V., Cohen,T. and Colijn,C. (2019) Beyond the SNP Threshold: Identifying Outbreak Clusters Using Inferred Transmissions. Mol. Biol. Evol., 36, 587–603.
