'''Script to randomly generate a metagenomic dataset.
Uses mason2 as subprocess.'''

import os
import subprocess
import argparse 
import ast 

import random
import math

from more_itertools import chunked

#################################### configuration ####################################
parser = argparse.ArgumentParser(description="Generate random dna sequences with singular genomes and parent genomes with children.",
                                 fromfile_prefix_chars='@')

parser.add_argument("output_dir", help="The directory where all output files are placed.")

parser.add_argument("-m", "--mason", default="", help="The path to the mason2 binary dir.")
parser.add_argument("-s", "--singulars", required=True, help="Sizes of genomes that are generated once.")
parser.add_argument("-p", "--parents", required=True, help="Sizes of genomes that are generated once and have children.")
parser.add_argument("-c", "--children", required=True, type=int, help="Number Of children per parent.")
parser.add_argument("-x", "--snp", default=0.001, type=float, help="Snp rate for generation of children.")
parser.add_argument("-i", "--indel", default=0.00001, type=float, help="Small indel rate for generation of children.")
parser.add_argument("-r", "--random-seeds", default=None,
                    help="Random seeds to use for all random processes (can be extracted from config_summary.txt of previous runs).")

args = parser.parse_args()

OUTPUT_DIR = args.output_dir
SINGULAR_FASTA_DIR = OUTPUT_DIR + "singular_genomes_fasta/"
PARENT_FASTA_DIR = OUTPUT_DIR + "parent_genomes_fasta/"
CHILD_VCF_DIR = OUTPUT_DIR + "child_genomes_vcf/"
CHILD_FASTA_DIR = OUTPUT_DIR + "child_genomes_fasta/"

MASON_DIR = args.mason

# genome generation
SINGULAR_GENOME_SIZES = ast.literal_eval(args.singulars.strip(' '))
PARENT_GENOME_SIZES = ast.literal_eval(args.parents.strip(' '))

# mason_variator
CHILD_GENOMES_PER_PARENT = args.children
SNP_RATE = args.snp
SMALL_INDEL_RATE = args.indel

# seed management
SEEDS_GIVEN = not args.random_seeds is None
RANDOM_SEEDS = ast.literal_eval(args.random_seeds.strip(' ')) if SEEDS_GIVEN else []
seed_index = 0

def next_random():
    global seed_index
    # generate new seeds if none were given
    if not SEEDS_GIVEN:
        RANDOM_SEEDS.append(int.from_bytes(os.urandom(3), byteorder="big"))
        return RANDOM_SEEDS[-1]
    
    # else use given seeds
    else:
        seed_index += 1
        return RANDOM_SEEDS[seed_index - 1]

random.seed(next_random())

#################################### execution ####################################

def check_error(proc, name):
    try:
        proc.check_returncode()
    except:
        print(f"---------- {name} failed with the following error output: ----------\n")
        print(proc.stderr.decode("utf-8"))
        quit()

# function for well padded filenames
total_genomes = len(SINGULAR_GENOME_SIZES) + (1 + CHILD_GENOMES_PER_PARENT) * len(PARENT_GENOME_SIZES)
bits = math.ceil(math.log10(total_genomes))
number_fmt = lambda i: f'{i:0{bits}d}'

# create all directories
os.mkdir(OUTPUT_DIR)
os.mkdir(SINGULAR_FASTA_DIR)
os.mkdir(PARENT_FASTA_DIR)
os.mkdir(CHILD_VCF_DIR)
os.mkdir(CHILD_FASTA_DIR)

def write_random_seq(filename, name, size):
    '''Write a random dna sequence with given size to a fasta file with given filename'''
    with open(filename, "w+") as f:
        f.write(">" + name + "\n")
        seq_gen = (random.choice(('A', 'C', 'G', 'T')) for _ in range(size))
        for seq_line in chunked(seq_gen, 80):
            f.write("".join(seq_line) + "\n")

fasta_file_listing = ""
completed_processes = []

# generate singular genomes
for i, size in enumerate(SINGULAR_GENOME_SIZES):
    name = "singular_" + number_fmt(i)
    filename = SINGULAR_FASTA_DIR + name + ".fasta"
    fasta_file_listing += filename + '\n'

    write_random_seq(filename, name, size)

parent_filenames = []
# generate parent genomes
for i, size in enumerate(PARENT_GENOME_SIZES):
    name = "init_" + number_fmt(i)
    filename = PARENT_FASTA_DIR + name + ".fasta"
    parent_filenames.append(filename)
    fasta_file_listing += filename + '\n'

    write_random_seq(filename, name, size)

# create child genomes with mason_variate
for parent, parent_filename in enumerate(parent_filenames):
    for child in range(CHILD_GENOMES_PER_PARENT):
        vcf_filename = CHILD_VCF_DIR + "child_" + number_fmt(parent) + "_" + number_fmt(child) + ".vcf"
        fasta_filename = CHILD_FASTA_DIR + "child_" + number_fmt(parent) + "_" + number_fmt(child) + ".fasta"

        fasta_file_listing += fasta_filename + '\n'

        proc = subprocess.run(
            [
            MASON_DIR + "mason_variator", "--verbose",
            "--seed", str(next_random()),
            "--snp-rate", str(SNP_RATE),
            "--small-indel-rate", str(SMALL_INDEL_RATE), 
            "-ir", parent_filename,
            "-ov", vcf_filename,
            "-of", fasta_filename,
            ],
            capture_output=True
        )

        check_error(proc, "mason_variate")
        completed_processes.append(proc)

# write files with mason output and error
stdout = "".join(map(lambda proc: proc.stdout.decode("ascii"), completed_processes))
stderr = "".join(map(lambda proc: proc.stderr.decode("ascii"), completed_processes))

with open(OUTPUT_DIR + "mason_stdout.txt", "w+") as f:
    f.write(stdout)

with open(OUTPUT_DIR + "mason_stderr.txt", "w+") as f:
    f.write(stderr)

# write file that lists all generated fasta files
with open(OUTPUT_DIR + "fasta_file_listing.txt", "w+") as f:
    f.write(fasta_file_listing)

# write info file with some info about what was generated
config_summary = (
    "------------------------------------------------------------------------\n"
    "Generated with the python script 'generate_dataset.py' which uses mason2\n"
    "\nConfigurations:\n\n"
    f"OUTPUT_DIR = {OUTPUT_DIR}\n"
    f"SINGULAR_GENOME_SIZES = {SINGULAR_GENOME_SIZES}\n"
    f"PARENT_GENOME_SIZES = {PARENT_GENOME_SIZES}\n"
    f"CHILD_GENOMES_PER_PARENT = {CHILD_GENOMES_PER_PARENT}\n"
    f"SNP_RATE = {SNP_RATE}\n"
    f"SMALL_INDEL_RATE = {SMALL_INDEL_RATE}\n"
    f"RANDOM_SEEDS = {RANDOM_SEEDS}\n"
    "------------------------------------------------------------------------\n"
)

with open(OUTPUT_DIR + "config_summary.txt", "w+") as f:
    f.write(config_summary)

print("Succesfully generated the dataset to " + OUTPUT_DIR)
