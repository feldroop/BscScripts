# BscScripts - Scripts for my bachelor thesis

The goal of this repository is to enable everyone to reproduce the tests and benchmarks I did for my bachelor thesis. All of the used scripts are written in Python. Feel free to adjust them if something does not work the way you like.

## Prerequisities

* Linux-like system. I use Ubuntu 20.04.
* Python 3. I use version 3.8.5.
* GCC. I use version 9.3.0.
* CMake. I use version 3.16.3.
* The Python library `more_itertools` (`pip install more_itertools`). Only needed for the `generate_datasets.py` script.
* The Python library `psutils` (`pip install psutils`) Only needed for the `compare.py` script.
* The SeqAn application `mason2`, download [here](http://packages.seqan.de/mason2/). Only needed for the `generate_datasets.py` script.
* Dependencies of `genome_updater`, see [here](https://github.com/pirovc/genome_updater). Only needed to download the real dataset.

## 1. Build Chopper binaries

First of all, clone [my fork of Chopper](https://github.com/Felix-Droop/Chopper) and build with `cmake`.

```
git clone --recurse-submodules https://github.com/Felix-Droop/Chopper
mkdir chopper_build
cd chopper_build
cmake ../Chopper
```

You can run the test to make sure everything works.
```
make test
```

## 2. Generate datasets

Clone this repository.

The `generate_dataset.py` script randomly generates a dataset of dna sequences. For more details see the help menu of the script. You can pass the parameters on the command line or write them into a file and pass the filename with `@`.

```
git clone https://github.com/Felix-Droop/BscScripts

python generate_dataset.py --help
python generate_dataset.py @config/generate_dataset_example.config
```

The script will generate a file `fasta_file_listing.txt` with all the names of sequence files in it. This will be needed later.

**Warning:** If you write the parameters into a file, make sure to place every single argument into a seperate line and have no trailing whitespaces. See examples.

## 3. Download real datasets

Clone the `genome_updater` repository. Run the tests script if you want to make sure it works.

```
git clone https://github.com/pirovc/genome_updater.git
cd genome_updater
source tests/tests.sh
```

Look at the README of the repo or at the help menu if you want to adjust the dataset. Otherwise you can use these parameters.

```
genome_updater.sh \
-d "refseq" \
-g "archaea,bacteria" \
-c "all" \
-l "Complete Genome" \
-f "genomic.fna.gz" \
-t 16 \
-o "/path/to/dataset/dir/"
```

**Warning:** This compressed dataset is roughly 20 GB large.

The script will not generate a sequence listing file. Therefore we must generate it ourselves.

```
python make_seq_list.py --help
python make_seq_list.py /path/to/files/ /output/file.txt
```

You can also specify the `--max-number` option to limit the number of files included. See help menu of the script for more details.

## 4. Run comparison

Now we can compare the different modes of `chopper pack` on our datasets. An example config file can be found in the config folder.

```
python compare.py --help
python compare.py @config/compare_example.config
```

The file will print a summary to the command line and also write it to a logfile in the output directory. All other outputs of subprocesses (`chopper pack`, etc.) are saved in that directory as well.
