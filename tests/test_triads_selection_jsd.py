import pytest
import numpy as np
import json
from cogent3.util.deserialise import deserialise_object
from clock_project.genome_analysis.traids_selection import extract_species_names, extract_taxonomic_order_info, get

from clock_project.replaced.triads_selection_jsd import (
    pseudorandom_ingroup_selection_between_order_jsd, 
    triads_binning_jsd_check, 
    get_outgroup_species)

@pytest.fixture
def sequence_path():
    return '/Users/gulugulu/repos/PuningAnalysis/results/output_data/homologies_sequence_common_name/ENSG00000162620.json'

@pytest.fixture
def sequences(sequence_path):
    with open (sequence_path, 'r') as infile:
        sequence_serialised = json.load(infile)
        sequence_test = deserialise_object(sequence_serialised)
        seqs_subset = sequence_test.take_seqs(['Eurasian_red_squirrel-ENSSVLG00005018133', 'Northern_American_deer_mouse-ENSPEMG00000000036', 'Vervet_AGM-ENSCSAG00000001381', 'Bolivian_squirrel_monkey-ENSSBOG00000023453', 'Dog-ENSCAFG00845004764', 'Bolivian_squirrel_monkey-ENSSBOG00000023453'])
    return seqs_subset

@pytest.fixture
def taxa_reference():
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_taxa_info.json', 'r') as infile:
        taxa_reference = json.load(infile)
    return taxa_reference

@pytest.fixture
def expected_species_names():
    return ['Eurasian_red_squirrel',
 'Northern_American_deer_mouse',
 'Vervet_AGM',
 'Bolivian_squirrel_monkey',
 'Dog',
 'Bolivian_squirrel_monkey']

@pytest.fixture
def expected_taxa_info():
    return {'Eurasian_red_squirrel': {'order': 'Rodentia',
   'family': 'Sciuridae',
   'genus': 'Sciurus'},
  'Northern_American_deer_mouse': {'order': 'Rodentia',
   'family': 'Cricetidae',
   'genus': 'Peromyscus'},
  'Vervet_AGM': {'order': 'Primates',
   'family': 'Cercopithecidae',
   'genus': 'Chlorocebus'},
  'Bolivian_squirrel_monkey': {'order': 'Primates',
   'family': 'Cebidae',
   'genus': 'Saimiri'},
  'Dog': {'order': 'Carnivora', 'family': 'Canidae', 'genus': 'Canis'}}

@pytest.fixture
def expected_order_info():
    return {'Rodentia': ['Eurasian_red_squirrel', 'Northern_American_deer_mouse'],
  'Primates': ['Vervet_AGM',
   'Bolivian_squirrel_monkey',
   'Bolivian_squirrel_monkey'],
  'Carnivora': ['Dog']}

@pytest.fixture
def order_reference():
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_order_info.json', 'r') as infile:
        order_reference = json.load(infile)
    return order_reference


def test_pseudorandom_ingroup_selection_between_order_jsd(sequences):
    available_seq = extract_species_names(sequences)
    ingroup_species_name = pseudorandom_ingroup_selection_between_order_jsd(available_seq, taxa_reference, order_reference, sequences)
    assert expected_taxa_info[ingroup_species_name['ingroup1']]['order'] != expected_taxa_info[ingroup_species_name['ingroup2']]['order']


