import random 
from more_itertools import chunked

def seed(a):
    random.seed(a)

def write_random_dna_seq_fasta(length, seq_id, filepath, mode):
    '''Write a random dna sequence with given size to a fasta file with given filename and file mode (a+ or w+).'''
    with open(filepath, mode) as f:
        f.write(">" + seq_id + "\n")
        seq_gen = (random.choice(('A', 'C', 'G', 'T')) for _ in range(length))
        for seq_line in chunked(seq_gen, 80):
            f.write("".join(seq_line) + "\n")