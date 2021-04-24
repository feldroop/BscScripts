import argparse
import pathlib
import subprocess

import dna_seq_util

from more_itertools import interleave

import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors

import numpy as np 

#################################### command line interface ####################################
parser = argparse.ArgumentParser(description="Measure quality of the chopper HyperLogLog implementation.",
                                 fromfile_prefix_chars='@')

io = parser.add_argument_group("Input/Output")
io.add_argument("-o", "--fasta-output", type=pathlib.Path, 
    help="If given, random sequences are generated, written to this file and used for the analysis.")
io.add_argument("-i", "--fasta-input", type=pathlib.Path, 
    help="If given, the sequences in this file are used for the evaluations.")
io.add_argument("-t", "--tsv-file", required=True, type=pathlib.Path, 
    help="If --fasta-output or --fasta-input was given, this file is generated. Else it is only read and the contents are evaluated")
io.add_argument("-c", "--chopper-bin", type=pathlib.Path, help="The binary directory of chopper.")

generation = parser.add_argument_group("Sequence generation")
generation.add_argument("-l", "--length", type=int, action="append",
    help="Adds a measurement with this sequence length. Can be specified multiple times.")
generation.add_argument("-n", "--number-seqs", type=int, 
    help="The number of sequences to use for each experiment.")
generation.add_argument("-s", "--seed", type=int, help="The random seed to use.")

experiments = parser.add_argument_group("HyperLogLog experiments")
experiments.add_argument("-b", "--bits", action="append",
    help="Adds an integer value which is tested as HyperLogLog bits parameter. Can be specified multiple times.")
experiments.add_argument("-k", "--kmer-size", default="20", 
    help="The size of the k-mers. Default is 20.")

args = parser.parse_args()

#################################### sequence generation ####################################

fasta_file = None

# sequences should be generated
if args.fasta_output:

    if not args.length or not args.number_seqs:
        print("Must specify --length and --number-segs to generate sequences.")
        quit()

    print("Generating sequences...")

    dna_seq_util.seed(args.seed)

    id = 0
    for length in args.length:
        for _ in range(args.number_seqs):
            dna_seq_util.write_random_dna_seq_fasta(length, f"seq{id}", args.fasta_output, "w+" if id == 0 else "a")
            id += 1
    
    fasta_file = args.fasta_output

#################################### hyperloglog experiments ####################################

if fasta_file is None:
    # hyperloglog experiments should be performed
    if args.fasta_input:
        fasta_file = args.fasta_input 

if fasta_file:
    if not args.bits or not args.chopper_bin:
        print("Must specify --bits and --chopper-bin to perform hyperloglog experiments.")
        quit()

    print("Building HyperLogLog sketches...")

    proc = subprocess.run(
        [
            args.chopper_bin / "measure_hyperloglog",
            "-i", fasta_file,
            "-o", args.tsv_file,
            "-k", args.kmer_size
        ] + list(interleave(["-b"] * len(args.bits), args.bits)),
        capture_output=True,
        encoding="utf-8"
    )

    if proc.returncode != 0:
        message = (
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}\n"
        )
        print(f"measure_hyperloglog failed with the following output:\n{message}")
        quit()

    else:
        print(f"measure_hyperloglog stdout:\n{proc.stdout}\n")
#################################### data analysis ####################################

print("Doing the evaluation...")

# expected fields are:
# sequence_id, sequence_length, sketch_register_size, estimated_cardinality,
# actual_cardinality, expected_relative_error, actual_relative_error

df = pd.read_csv(
    args.tsv_file,
    sep="\t",
    comment='#',
    header=0,
)

fig = plt.figure()
ax = plt.gca()

colors = {}
# add horizontal dotted lines for each register size for the control
for reg_size, group in df.groupby('sketch_register_size'):
    
    # take some color out of the tableau colors and fix it for the given register size
    colors[reg_size] = list(mcolors.TABLEAU_COLORS.keys())[len(colors)]

    plt.axhline(
        y=group['expected_relative_error'].iloc[0], 
        linestyle='dotted',
        label=f"{reg_size} registers/bytes",
        color=colors[reg_size]
    )

ax.set_ylim(0.0, 0.1)
ax.set_ylabel("Relative error")
ax.set_xlabel("Sequence length")

pos = 1
last_seq_length = None
for (seq_length, reg_size), group in df.groupby(['sequence_length','sketch_register_size']):

    # if seq_length changed, we want one empty position for visual seperation
    if last_seq_length and last_seq_length != seq_length:
        pos += 1
    last_seq_length = seq_length

    boxplot = ax.boxplot(
        group['actual_relative_error'],
        positions=[pos],
        labels=['{:.1e}'.format(seq_length)],
        widths=[0.5]
    )

    # color the boxplots
    for lines in boxplot.values():
        for line in lines:
            line.set_color(colors[reg_size])
            line.set_linewidth(1.5)
        
    pos += 1
    
plt.legend()
plt.show()
