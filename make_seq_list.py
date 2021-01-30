import os 
import argparse
import math
import random 

DATA_DIR = "/home/felix/data/bsc_thesis/real_archaea_dataset/2021-01-27_20-12-17/files/"
OUT_FILE = "/home/felix/data/bsc_thesis/real_archaea_dataset/fasta_file_listing.txt"

parser = argparse.ArgumentParser(description="Generate a file with the names of all files inside a folder. Needed for chopper pack.",
                                 fromfile_prefix_chars='@')

parser.add_argument("data_dir", help="The directory where the files are stored.")
parser.add_argument("out_file", help="The location where the output file should be created")
parser.add_argument("-m", "--max-number", type=int, default=math.inf, 
                    help="The maximum number of files to include in the list. If not all are included, the subset is picked at random.")

args = parser.parse_args()

DATA_DIR = args.data_dir
OUT_FILE = args.out_file 
MAX_NUM = args.max_number

file_list = os.listdir(DATA_DIR)

if len(file_list) <= MAX_NUM:
    full_paths = map(lambda f: DATA_DIR + f, file_list)
    seq_listing = "\n".join(full_paths)

else:
    random.shuffle(file_list)

    seq_listing = ""
    for filename in file_list[:MAX_NUM]:
        seq_listing += DATA_DIR + filename + "\n"

with open(OUT_FILE, "w+") as f:
    f.write(seq_listing)
