# PuningAnalysis
## Data sampling 

Sequences from Ensembl release 113 raw_data

### Eliminate low-quality sequences 
based on sequence length, matching proportion to human reference, and proportion of degenerate position.
data_filter_seqs.py
input: raw_data -> output: sampled_unaligned

### Aligning the sequences
amino acid translation MAFFT
