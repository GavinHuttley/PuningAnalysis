from cogent3 import open_
import os
import glob
import pathlib
from collections import Counter
from cogent3.parse.fasta import MinimalFastaParser

#directory that contain all species folders
base_dir = "/home2/puning/honours/data/whole_genome_mammal87/ensembl_download/genomes"

all_fasta_paths = {}

for species in os.listdir(base_dir):
    fasta_dir = os.path.join(base_dir, species, 'fasta')
    fasta_files = glob.glob(os.path.join(fasta_dir, '*.fa.gz'))
    all_fasta_paths[species] = fasta_files


def counts_nucs(path):
    with open_(path) as infile:
        data = infile.readlines()

    parser = MinimalFastaParser(data)
    nuc_counts = Counter()
    for label, seq in parser:
        nuc_counts.update(seq)
        del seq # delete the sequence instance

    # discard non-canonical nucleotides
    return {b: nuc_counts[b] for b in "ACGT"}

nuc_counts = {}
for species, paths in all_fasta_paths.items():
    print(species)
    nuc_counts[species] = [counts_nucs(path) for path in paths]


print(nuc_counts)
