import argparse
import subprocess
import pathlib 
import os 
import time 

import collections
from enum import Enum 

import pandas as pd

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
parser.add_argument("-q", "--quick", action="store_true", help="If given, assume that the binning file already exists.")

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

if not args.quick:
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


    print_and_log(
        f"---------- multilevel packing done. ----------\n"
        f"           took {round(elapsed_time, 3)} seconds.\n"
    )
else:
    print_and_log("---------- skipped execution. ----------\n")

#################################### execution ####################################
df = pd.read_csv(
    binning_filename, 
    sep="\t", 
    comment="#", 
    header=None,
    names=["File", "Bin_Index", "Num_Bins", "Cardinality_Estimate"]
)

# represents a technical bin
class Bin:
    class Type(Enum):
        Split = 0
        Merged = 1
    
    def __init__(self):
        self.type = self.Type.Split
        self.cardinality_estimate = 0
        self.contained_ubs = 0
        self.num_bins = 0
        self.child_bins = collections.defaultdict(lambda: Bin())

top_level_bins = collections.defaultdict(lambda: Bin())

def to_tup(s):
    return tuple(map(int, s.split(";")))

# extract data from flat df and turn into useful hierarchical structure
for _, (_, *info) in df.iterrows():
    bin_index, num_bins, cardinality_estimate = tuple(map(to_tup, info))

    curr_level_bins = top_level_bins
    for level in range(len(bin_index)):
        curr_bin = curr_level_bins[bin_index[level]]

        curr_bin.contained_ubs += 1
        if curr_bin.contained_ubs > 1:
            curr_bin.type = Bin.Type.Merged
        
        curr_bin.cardinality_estimate = cardinality_estimate[level]
        curr_bin.num_bins = num_bins[level]
        
        curr_level_bins = curr_bin.child_bins

# statistics for a level of the HIBF
class Statistics():
    def __init__(self):
        self.num_ibs = 0
        self.num_bins = 0
        self.split_bins = 0
        self.merged_bins = 0
        self.num_split_ubs = 0
        self.num_merged_ubs = 0
        self.max_ubs_in_split = 0
        self.max_ubs_in_merged = 0
        self.s_tech = 0
    
    def __str__(self):
        return (
            f"\n"
            f"#IBFS:                 : {self.num_ibs:,}\n"
            f"#bins                  : {self.num_bins:,}\n"
            f"#split bins            : {self.split_bins:,}\n"
            f"#merged bins           : {self.merged_bins:,}\n"
            f"#UBs in split bins     : {self.num_split_ubs:,}\n"
            f"#UBs in merged bins    : {self.num_merged_ubs:,}\n"
            f"max #UBs in split bin  : {self.max_ubs_in_split:,}\n"
            f"max #UBs in merged bin : {self.max_ubs_in_merged:,}\n"
            f"Total S_tech           : {self.s_tech:,}\n"
        )

levels = collections.defaultdict(lambda: Statistics())

# recursively gather statistics for a given set of bins and a given level
def gather_statistics(level, bins):
    stat = levels[level]
    stat.num_ibs += 1

    max_bin_card = 0
    local_num_bins = 0

    for bin in bins.values():
        stat.num_bins += bin.num_bins
        local_num_bins += bin.num_bins

        if bin.type == Bin.Type.Split:
            stat.split_bins += bin.num_bins

            stat.num_split_ubs += bin.contained_ubs
            stat.max_ubs_in_split = max(stat.max_ubs_in_split, bin.num_bins)
        
        else:
            stat.merged_bins += bin.num_bins
            stat.num_merged_ubs += bin.contained_ubs
            stat.max_ubs_in_merged = max(stat.max_ubs_in_merged, len(bin.child_bins))
            gather_statistics(level + 1, bin.child_bins)
        
        max_bin_card = max(max_bin_card, bin.cardinality_estimate)

    stat.s_tech += max_bin_card * local_num_bins

# gather all statistics
gather_statistics(0, top_level_bins)

# print and log statistics for all levels
for level, stat in levels.items():
    print_and_log(f"Level {level}:\n{stat}")

print_and_log("Total S_tech = sum over all IBFs on the given level of (#bins * <maximum bin cardinality>)\n")
