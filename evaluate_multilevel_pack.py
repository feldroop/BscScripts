import argparse
import subprocess
import pathlib 
import os 
import time 

import pandas as pd
import numpy as np

# timestamp
t = time.localtime()
timestamp = f"{t.tm_year}-{t.tm_mon}-{t.tm_mday}_{t.tm_hour}-{t.tm_min}-{ t.tm_sec}"

#################################### configuration ####################################
parser = argparse.ArgumentParser(description="Run the multilevel pack algorithm and evaluate the results",
                                 fromfile_prefix_chars='@')

parser.add_argument("-o", "--output-dir", required=True, type=pathlib.Path,
                    help="The directory where all output files are placed.")
parser.add_argument("-k", "--kmer-count-file", required=True, type=pathlib.Path, 
                    help="The file for chopper pack in which all sequence files are listed with k-mer counts.")
parser.add_argument("-c", "--chopper-bin-dir", required=True, type=pathlib.Path, help="The binary directory of chopper.")
parser.add_argument("-d", "--hll-dir", required=True, help="The dir where the hlls are cached.", type=pathlib.Path)
parser.add_argument("-l", "--log", default=f"{timestamp}_log.txt", help="The name for the log file (not the whole path).")
parser.add_argument("-b", "--bins", required=True, type=int, help="The number of technical bins for chopper pack.")
parser.add_argument("-a", "--alpha", default=1.2, type=float, help="The alpha for the internal binning algorithm.")
parser.add_argument("-m", "--max-ratio", default=0.0, type=float, 
                    help="The maximal cardinality ratio in the clustering intervals (must be < 1).")
parser.add_argument("-t", "--threads", default=1, type=int, help="The number of threads to use.")
parser.add_argument("-p", "--false-positive-rate", default=0.05, type=float, help="The false positive rate for the IBFs.")
parser.add_argument("-s", "--num-hash-functions", default=2, type=int, help="The number hash functions for the IBFs.")

args = parser.parse_args()

#################################### execution ####################################
if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)

# setup logging
log_path = args.output_dir / args.log

def print_and_log(message):
    print(message)
    with open(log_path, "a+") as f:
        f.write(message + '\n')

print_and_log(
    "\n---------- configuration: ----------\n\n"
    f"output directory: {args.output_dir}\n"
    f"k-mer counts    : {args.kmer_count_file}\n"
    f"chopper binaries: {args.chopper_bin_dir}\n"
    f"hll directory   : {args.hll_dir}\n"
    f"log file        : {log_path}\n\n"
    f"pack bins       : {args.bins}\n"
    f"pack alpha      : {args.alpha}\n"
    f"max ratio       : {args.max_ratio}\n"
    f"threads         : {args.threads}\n"
    f"FPR             : {args.false_positive_rate}\n"
    f"Hash functions  : {args.num_hash_functions}\n"
)

binning_filename = args.output_dir / "multilevel.binning"
output_filename = args.output_dir / "pack_multilevel_full_output.txt"

start_time = time.perf_counter()

pack_proc = subprocess.run([
    args.chopper_bin_dir / "chopper", 
    "pack",
    "-f", args.kmer_count_file,
    "-d", args.hll_dir,
    "-b", str(args.bins),
    "-a", str(args.alpha),
    "-m", str(args.max_ratio),
    "-t", str(args.threads),
    "-p", str(args.false_positive_rate),
    "-s", str(args.num_hash_functions),
    "-o", binning_filename
    ],
    encoding='utf-8',
    capture_output=True
)

elapsed_time = time.perf_counter() - start_time

message = (
        f"---------- stdout ----------\n"
        f"{pack_proc.stdout}\n"
        f"---------- stderr ----------\n"
        f"{pack_proc.stderr}\n"
)

if pack_proc.returncode != 0:
    print_and_log(f"---------- multilevel pack failed with the following output: ----------\n\n"
                  f"{message}"
    )
    quit()

with open(output_filename, "w+") as f:
    f.write(message)

for line in pack_proc.stdout.splitlines():
    if 'optimum' in line:
        print_and_log(
            f"---------- multilevel packing done. {line} ----------\n"
            f"           took {round(elapsed_time, 3)} seconds."
        )

for line in pack_proc.stderr.splitlines():
    if 'peak memory usage' in line:
        print_and_log("           " + line + '\n')

#################################### execution ####################################
df = pd.read_csv(
    binning_filename, 
    sep="\t", 
    comment="#", 
    header=None,
    names=["File", "Bin_Index", "Num_Bins", "Est_Max_TB"]
)

print(df)
