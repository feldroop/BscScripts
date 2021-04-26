import os
import argparse
import pathlib 
import time
import subprocess

# timestamp
t = time.localtime()
timestamp = f"{t.tm_year}-{t.tm_mon}-{t.tm_mday}_{t.tm_hour}-{t.tm_min}-{ t.tm_sec}"

#################################### configuration ####################################
parser = argparse.ArgumentParser(description="Count k-mers from files, run different packing algorithms and evaluate results",
                                 fromfile_prefix_chars='@')

parser.add_argument("output_dir", help="The directory where all output files are placed.", 
                    type=pathlib.Path)
parser.add_argument("seqfile_list_file", help="The file for chopper pack in which all sequence files are listed.",
                    type=pathlib.Path)
parser.add_argument("binary_dir", help="The binary directory of chopper.", 
                    type=pathlib.Path)

parser.add_argument("-d", "--hll-dir", required=True, help="The dir where the hlls are cached.",
                    type=pathlib.Path)
parser.add_argument("-l", "--log", default=f"{timestamp}_log.txt", help="The name for the log file (not the whole path).")
parser.add_argument("-b", "--bins", required=True, type=int, help="The number of technical bins for chopper pack.")
parser.add_argument("-k", "--kmer-size", default=20, type=int, help="The size of the k-mers.")
parser.add_argument("-a", "--alpha", default=1.2, type=float, help="The alpha for the internal binning algorithm.")
parser.add_argument("-s", "--sketch-bits", default=12, type=int, 
                    help="The number of bits to distribute values for the HyperLogLog sketches.")
parser.add_argument("-m", "--max-ratio", default=0.5, type=float, 
                    help="The maximal cardinality ratio in the clustering intervals (must be < 1).")
parser.add_argument("-t", "--threads", default=1, type=int, help="The number of threads to use.")
parser.add_argument("-x", "--no-recount", action='store_true', 
                    help="If given, chopper count is not invoked and kmer_counts.txt from output dir is used.")
parser.add_argument("-p", "--peak-memory", action='store_true',
                    help="If given, the peak memory usage of chopper processes is computed. Might be a slight slow-down.")
parser.add_argument("-e", "--exclusively-hlls", action='store_true',
                    help="If given, the hll counts are used for chopper pack instead of the eact counts.")

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
    f"seqfile list    : {args.seqfile_list_file}\n"
    f"chopper binaries: {args.binary_dir}\n"
    f"hll directory   : {args.hll_dir}\n"
    f"log file        : {log_path}\n\n"
    f"pack bins  : {args.bins}\n"
    f"k-mer size : {args.kmer_size}\n"
    f"pack alpha : {args.alpha}\n"
    f"sketch bits: {args.sketch_bits}\n"
    f"max ratio  : {args.max_ratio}\n"
    f"threads    : {args.threads}\n"
    f"no recount : {args.no_recount}\n"
    f"hll counts : {args.exclusively_hlls}\n"
)

def handle_outputs(proc, name, filename):
    '''If the process errored, print the error, else write stdout to a file'''

    message = (
        f"---------- stdout ----------\n"
        f"{proc.stdout}\n"
        f"---------- stderr ----------\n"
        f"{proc.stderr}\n"
    )

    if proc.returncode != 0:
        print_and_log(f"---------- {name} failed with the following output: ----------\n\n"
                      f"{message}"
        )
        quit()
    
    with open(filename, "w+") as f:
        f.write(message)

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

def run_count(extra_flags, name):
    kmer_counts_filename = args.output_dir / (name + "_kmer_counts.txt")
    
    start_time = time.perf_counter()

    count_proc = subprocess.run([
        args.binary_dir / "chopper", 
        "count",
        "-f", args.seqfile_list_file,
        "-o", kmer_counts_filename,
        "-k", str(args.kmer_size),
        "-t", str(args.threads),
        "-s", str(args.sketch_bits),
        "--disable-minimizers",
        ] + extra_flags,
        encoding='utf-8',
        capture_output=True
    )
    
    elapsed_time = time.perf_counter() - start_time

    output_filename = args.output_dir / (name + "_count_outputs.txt")
    handle_outputs(count_proc, "chopper count", output_filename)

    print_and_log(
        f"---------- k-mer counting with {name} counts done ----------\n"
        f"           took {round(elapsed_time, 3)} seconds."
    )

    for line in count_proc.stdout.splitlines():
        if 'peak memory usage' in line:
            print_and_log("           " + line + '\n')

def run_pack(extra_flags, name):
    kmer_counts_filename = args.output_dir / (("hll" if args.exclusively_hlls else "exact") + "_kmer_counts.txt")
    binning_filename = args.output_dir / (name + ".binning")
    output_filename = args.output_dir / ("pack_" + name + "_full_output.txt")

    start_time = time.perf_counter()

    pack_proc = subprocess.run([
        args.binary_dir / "chopper", 
        "pack",
        "-f", kmer_counts_filename,
        "-d", args.hll_dir,
        "-b", str(args.bins),
        "-a", str(args.alpha),
        "-m", str(args.max_ratio),
        "-o", binning_filename
        ] + extra_flags,
        encoding='utf-8',
        capture_output=True
        )
    
    elapsed_time = time.perf_counter() - start_time

    handle_outputs(pack_proc, f"chopper pack with {name}", output_filename)

    for line in pack_proc.stdout.splitlines():
        if 'optimum' in line:
            print_and_log(
                f"---------- packing with {name} done. {line} ----------\n"
                f"           took {round(elapsed_time, 3)} seconds."
            )
        elif 'peak memory usage' in line:
            print_and_log("           " + line + '\n')

def evaluate(name):
    kmer_counts_filename = args.output_dir / "exact_kmer_counts.txt"
    binning_filename = args.output_dir / (name + ".binning")

    proc = subprocess.run([
        args.binary_dir / "count_HIBF_kmers_based_on_binning", 
        "-f", binning_filename,
        "-c", kmer_counts_filename,
        "-k", str(args.kmer_size),
        "-t", str(args.threads)
        ],
        encoding='utf-8',
        capture_output=True
    )

    output_filename = args.output_dir / f"evaluation_{name}.txt"
    handle_outputs(proc, f"count_HIBF_kmers_based_on_binning for the {name}", output_filename)

    maxi, splits, merges, low_level_size = analyze_result(proc.stdout)
    print_and_log(
        f"---------- evaluating with {name} done. ----------\n\n"
        f"#split bins              : {splits}\n"
        f"#merged bins             : {merges}\n"
        f"largest bin              : {maxi}\n"
        f"largest bin * #bins      : {maxi * args.bins}\n"
        f"lower level k-mers (sum) : {low_level_size}\n"
        f"\n{proc.stdout if len(proc.stdout.splitlines()) <= 64 else ''}"
        )

if not args.no_recount:
    # run chopper count on the fasta listing
    run_count([], "exact")
    run_count(["-e", "-d", str(args.hll_dir)], "hll")

else:
    print_and_log("---------- No recount of k-mers done. ----------\n")

# run chopper pack WITHOUT union estimates
run_pack([], "reference")

# run chopper pack WITH union estimates
run_pack(["-u"], "union")

# run chopper pack WITH union estimates AND rearranging
run_pack(["-u", "-r"], "rearrange")

# run count_HIBF_kmers_based_on_binning for the reference, unions and rearrange result
evaluate("reference"), 
evaluate("union"), 
evaluate("rearrange")
