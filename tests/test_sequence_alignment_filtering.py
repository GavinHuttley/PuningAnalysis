import pytest
import numpy as np
from clock_project.genome_analysis.sequence_alignment_filtering import drop_low_matching, get_score_at_quantile, get_matching_score_to_human, drop_ambiguous, seq_without_redundent, common_name_renamer, drop_short_seq
from cogent3 import get_app, open_data_store

loader = get_app("load_unaligned", format="fasta", moltype="dna")

@pytest.fixture
def seqs_path():
    seq_dir = "/Users/gulugulu/Desktop/honours/data_local_2/sampled_homologies"
    input_dstore_seqs = open_data_store(seq_dir, suffix="fa", mode="r")
    path = input_dstore_seqs.completed[0]
    return path


@pytest.fixture
def seqs(seqs_path):
    return loader(seqs_path)

@pytest.fixture
def sub_seqs(seqs):
    return seqs.take_seqs(['echinops_telfairi-ENSETEG00000000236', 'dasypus_novemcinctus-ENSDNOG00000018355', 'loxodonta_africana-ENSLAFG00000013661'])

@pytest.fixture
def common_names():
    return ['Lesser_hedgehog_tenrec', 'Armadillo', 'Elephant']

@pytest.fixture
def quantile1():
    return 0.05

@pytest.fixture
def quantile2():
    return 0.02

def test_drop_low_matching(seqs, quantile1):
    seq_scores = get_matching_score_to_human(seqs)
    score_at_quantile = get_score_at_quantile(seq_scores, quantile1)
    seqs_filtered = drop_low_matching(seqs)
    seq_scores_after = get_matching_score_to_human(seqs_filtered)
    for score in seq_scores_after.values():
        assert score > score_at_quantile

def test_seq_without_redundent(seqs):
    seqs_filtered = seq_without_redundent(seqs)
    seq_names = seqs_filtered.names
    genus_list = []
    for name in seq_names:
        genus = name.split('_')[0]
        genus_list.append(genus)
    assert len(set(genus_list)) == len(genus_list)

def test_drop_ambiguous(seqs, quantile2):
    seqs_filtered = drop_ambiguous(seqs)
    seqs_filtered_lengths_wa = seqs_filtered.get_lengths(include_ambiguity=False)
    seqs_filtered_lengths_woa = seqs_filtered.get_lengths(include_ambiguity=True)
    for name in seqs_filtered_lengths_wa.keys():
        assert seqs_filtered_lengths_woa[name] / seqs_filtered_lengths_wa[name] >= 1 - quantile2


def test_common_name_renamer(sub_seqs, common_names):
    sub_seqs_renamed = common_name_renamer(sub_seqs)
    print(sub_seqs_renamed)
    seq_names = sub_seqs_renamed.names
    for name in seq_names:
        assert name in common_names
    for name in common_names:
        assert name in seq_names

def test_drop_short_seq(seqs):
    seqs_filtered = drop_short_seq(seqs)
    for species, length in seqs_filtered.get_lengths().items():
        assert length >= 500

