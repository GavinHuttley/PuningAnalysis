# PuningAnalysis
## Data sampling 

Sequences from Ensembl release 113 raw_data

### Eliminate low-quality sequences 
based on sequence length, matching proportion to human reference, and proportion of degenerate positions.
data_filter_seqs.py
input: raw_data -> output: sampled_unaligned

### Aligning the sequences
Amino acid translation MAFFT
mafft_align.py
input: sampled_unaligned -> output: raw_aligned

### Filtering aligned sequence
Filter the alignment by only keeping the third codon position, removing degenerate nucleotides, remove gene wiht less than 30 species and sequences with a  length of less than 550 nucleotides.
data_filter_aligns.py
input: raw_aligned -> sampled_aligned

### Fit the alignments to the GN model
Fitting the alignment for each gene to a general nucleotide model and getting a phylogenetic tree for each gene.
get_model_fitting_tree.py
input: sampled_aligned -> output: whole_gene_model_fitting

### Taxanonic triple selection
sample taxanomic triple for each gene based on the phylogenetic relation in the tree generated from the  last step, the raw sequence of the three species was aligned again using the same alignment process in MAFFT_align.py
triple_selection.py
input: whole_gene_model_fitting, sampled_unaligned, sampled_aligned -> taxanomic_triples, taxanomic_triples_alignments

### Triple model fitting and Taxanomic triple information collection
Fit the taxanomic triple alignment to the general nucleotide model and extract associated information for each taxanomic triple, including:
{'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'ingroup_jsd': ingroup_jsd,
        'nuc_freqs_dict': distributions,
        'nabla_values': nabla_dict,
        'matrices': matrices_list}

The triple information was saved into a JSON file ('triples_info_dict_new.json') for each gene in the associated folder in the triple model filter folder.
mammal_orthologs_hsap_1/triples_model_fitting/'gene_id'/triples_info_dict_new.json
mammal_orthologs_hsap_1/triples_model_fitting/'gene_id'/model_fitting_result

triple_model_fitting.py, triples_info_collection.py 
Input: taxanomic_triples_alignments -> ouput: triple_model_fitting.


### Taxanomic triple representative subset collection.
Stratified sampling taxanomic triple representative subset that contains the same proportion of triples from each gene compared to the full set. 
triples_representative_subset_selection.py
Input: triple_model_fitting, taxanomic_triples_alignments -> Oputput: triples_representative_subset, triples_representative_subset_json(selected triple information)

### Select valid substitution matrices
Collect matrices in the taxanomic triples where the rate parameter does not approach the boundary value. 
valid_matrix_collection.py, yapeng_check_BV.py
input: triples_model_fitting --> valid_matrix_full.json

### TOC hypothesis test result
toc_hypothesis_test_results

### TOS bootstrap result 
tos_boostrapping_results
