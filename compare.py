import os
import subprocess
import time
import argparse

#################################### configuration ####################################
parser = argparse.ArgumentParser(description="Count k-mers from files, run different packing algorithms and evaluate results",
                                 fromfile_prefix_chars='@')

parser.add_argument("output_dir", help="The directory where all output files are placed.")
parser.add_argument("seqfile_list_file", help="The file for chopper pack in which all sequence files are listed.")
parser.add_argument("binary_dir", help="The bianry directory of chopper.")

parser.add_argument("-c", "--hll-cache-dir", required=True, help="The dir where the hlls are cached.")
parser.add_argument("-b", "--bins", type=int, required=True, help="The number of technical bins for chopper pack.")
parser.add_argument("-k", "--kmers", default=20, type=int, help="The size of the k-mers.")
parser.add_argument("-a", "--alpha", default=1, type=int, help="The alpha for the internal binning algorithm.")
parser.add_argument("-s", "--sketch-bits", default=12, type=int, 
                    help="The number of bits to distribute values for the HyperLogLog sketches.")
parser.add_argument("-t", "--threads", default=1, type=int, help="The number of threads to use.")
parser.add_argument("-x", "--no-recount", action='store_true', 
                    help="If given, chopper count is not invoked and kmer_counts.txt from output dir is used.")

args = parser.parse_args()

OUTPUT_DIR = args.output_dir
SEQ_LIST = args.seqfile_list_file
BINARY_DIR = args.binary_dir
HLL_CACHE_DIR = args.hll_cache_dir

KMER_SIZE = args.kmers

PACK_BINS = args.bins
THREADS = args.threads
PACK_ALPHA = args.alpha
SKETCH_BITS = args.sketch_bits

NO_RECOUNT = args.no_recount
#################################### execution ####################################

if not os.path.isdir(OUTPUT_DIR):
    os.mkdir(OUTPUT_DIR)

# setup logging
t = time.localtime()
t_str = "-".join(map(str, (t.tm_year, t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec)))
log_filename = OUTPUT_DIR + t_str + "_log.txt"

def print_and_log(string):
    print(string)
    with open(log_filename, "a+") as f:
        f.write(string + '\n')

print_and_log(
    "---------- configuration: ----------\n\n"
    f"k-mer size : {KMER_SIZE}\n"
    f"pack bins  : {PACK_BINS}\n"
    f"pack alpha : {PACK_ALPHA}\n"
    f"sketch bits: {SKETCH_BITS}\n"
    f"threads    : {THREADS}\n"
)

def handle_outputs(proc, name, filename, with_stderr):
    '''If the process errored, print the error, else write stdout to a file'''
    if proc.returncode != 0:
        print_and_log(f"---------- {name} failed with the following output: ----------\n\n"
              f"---------- stdout ----------\n"
              f"{proc.stdout.decode('utf-8')}\n"
              f"---------- stderr ----------\n"
              f"{proc.stderr.decode('utf-8')}\n"
        )
        quit()
    
    output_str = proc.stdout.decode('utf-8')
    
    if with_stderr:
        output_str += proc.stderr.decode('utf-8')
    
    with open(filename, "w+") as f:
        f.write(output_str)

def analyze_result(s):
    '''Find the biggest technical bin from count_HIBF_kmers_based_on_binning output'''
    maxi, splits, merges, low_level_size = 0, 0, 0, 0

    for line in s.splitlines():
        split_line = line.split('\t')

        if len(split_line) > 1:
            maxi = max(maxi, int(split_line[1]))

        if "SPLIT_BIN" in split_line[0]:
            splits += 1

        if "MERGED_BIN" in split_line[0]:
            merges += 1
            low_level_size += int(split_line[2])    

    return (maxi, splits, merges, low_level_size)

kmer_counts_filename = OUTPUT_DIR + "kmer_counts.txt"

def run_pack(extra_flags, name):
    binning_filename = OUTPUT_DIR + name + ".binning"
    output_filename = OUTPUT_DIR + "pack_" + name + "_full_output.txt"

    pack_proc = subprocess.run([
        BINARY_DIR + "chopper", 
        "pack",
        "-f", kmer_counts_filename,
        "-c", str(HLL_CACHE_DIR),
        "-b", str(PACK_BINS),
        "-k", str(KMER_SIZE),
        "-a", str(PACK_ALPHA),
        "-s", str(SKETCH_BITS),
        "-t", str(THREADS),
        "-o", binning_filename
        ] + extra_flags,
        capture_output=True
        )

    handle_outputs(pack_proc, f"chopper pack with {name}", output_filename, True)

    for line in pack_proc.stdout.decode("utf-8").splitlines():
        if 'optimum' in line:
            print_and_log(f"---------- packing with {name} done. {line} ----------")

def evaluate(name):
    binning_filename = OUTPUT_DIR + name + ".binning"

    proc = subprocess.run([
        BINARY_DIR + "count_HIBF_kmers_based_on_binning", 
        "-f", binning_filename,
        "-k", str(KMER_SIZE),
        "-t", str(THREADS)
        ],
        capture_output=True
    )

    output_filename = OUTPUT_DIR + f"evaluation_{name}.txt"
    handle_outputs(proc, f"count_HIBF_kmers_based_on_binning for the {name}", output_filename, True)

    evaluation_output = proc.stdout.decode('utf-8')
    maxi, splits, merges, low_level_size = analyze_result(evaluation_output)
    print_and_log(
        f"---------- evaluating with {name} done. ----------\n\n"
        f"Number of split  bins: {splits}\n"
        f"Number of merged bins: {merges}\n"
        f"Maximum technical bin: {maxi}\n"
        f"High level k-mers    : {maxi * PACK_BINS}\n"
        f"Low  level k-mers    : {low_level_size}\n"
        f"Total k-mers (alpha) : {maxi * PACK_BINS + low_level_size * PACK_ALPHA}\n\n"
        f"{evaluation_output if len(evaluation_output.splitlines()) <= 64 else ''}"
        )

if not NO_RECOUNT:
    # run chopper count on the fasta listing
    count_proc = subprocess.run([
        BINARY_DIR + "chopper", 
        "count",
        "-f", SEQ_LIST,
        "-k", str(KMER_SIZE),
        "-t", str(THREADS),
        "--disable-minimizers"
        ],
        capture_output=True
        )
    
    handle_outputs(count_proc, "chopper count", kmer_counts_filename, False)

print_and_log("---------- k-mer counting done ----------\n")

# run chopper pack WITHOUT union estimates
run_pack([], "reference")

# run chopper pack WITH union estimates
run_pack(["-u"], "union")

# run chopper pack WITH union estimates AND resorting
run_pack(["-u", "-r"], "resort")

print_and_log('\n')

# run count_HIBF_kmers_based_on_binning for the reference, unions and resort result
evaluate("reference"), 
evaluate("union"), 
evaluate("resort")
