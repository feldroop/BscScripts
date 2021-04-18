import os 
import argparse
import random 
import pathlib 

parser = argparse.ArgumentParser(description="Generate a file with the names of all files inside a folder. Needed for chopper pack.",
                                 fromfile_prefix_chars='@')

parser.add_argument("data_dir", help="The directory where the files are stored.", type=pathlib.Path)
parser.add_argument("out_file", help="The location where the output file should be created", type=pathlib.Path)
parser.add_argument("-m", "--max-number", type=int, default=1000000000000, 
                    help="The maximum number of files to include in the list. If not all are included, the subset is picked at random.")

args = parser.parse_args()

file_list = os.listdir(args.data_dir)
random.shuffle(file_list)

if len(file_list) <= args.max_number:
    full_paths = map(lambda f: str(args.data_dir / f), file_list)
    seq_listing = "\n".join(full_paths)

else:
    seq_listing = ""
    for filename in file_list[:args.max_number]:
        seq_listing += str(args.data_dir / filename)  + "\n"

with open(args.out_file , "w+") as f:
    f.write(seq_listing)
