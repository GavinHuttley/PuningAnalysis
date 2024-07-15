import pytest
import numpy as np
import json
from cogent3.util.deserialise import deserialise_object
from clock_project.genome_analysis.traids_selection import extract_species_names, extract_taxonomic_order_info, pseudorandom_sequence_selection_within_order, pseudorandom_sequence_selection_between_order, initialize_bins, traids_selecting_binning, get_traits_info

@pytest.fixture
def sequence_path():
    return '/Users/gulugulu/repos/PuningAnalysis/results/output_data/homologies_sequence_common_name/ENSG00000117000.json'

@pytest.fixture
def sequences(sequence_path):
    with open (sequence_path, 'r') as infile:
        sequence_serialised = json.load(infile)
        sequence_test = deserialise_object(sequence_serialised)
        seqs_subset = sequence_test.take_seqs(['Gibbon-ENSNLEG00000012124', 'Hedgehog-ENSEEUG00000001288', 'Armadillo-ENSDNOG00000035279', 'Blue_whale-ENSBMSG00010003196', 'Olive_baboon-ENSPANG00000018124', 'Megabat-ENSPVAG00000015936'])
    return seqs_subset

@pytest.fixture
def taxa_reference():
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_taxa_info.json', 'r') as infile:
        taxa_reference = json.load(infile)
    return taxa_reference

@pytest.fixture
def expected_species_names():
    return ['Gibbon', 'Hedgehog', 'Armadillo', 'Blue_whale', 'Olive_baboon', 'Megabat']

@pytest.fixture
def expected_taxa_info():
    return {'Gibbon': {'order': 'Primates', 'family': 'Hylobatidae', 'genus': 'Various'},
            'Hedgehog': {'order': 'Eulipotyphla',
            'family': 'Erinaceidae',
            'genus': 'Erinaceus'},
            'Armadillo': {'order': 'Xenarthra',
            'family': 'Dasypodidae',
            'genus': 'Dasypus'},
            'Blue_whale': {'order': 'Cetacea',
            'family': 'Balaenopteridae',
            'genus': 'Balaenoptera'},
            'Olive_baboon': {'order': 'Primates',
            'family': 'Cercopithecidae',
            'genus': 'Papio'},
            'Megabat': {'order': 'Chiroptera',
            'family': 'Pteropodidae',
            'genus': 'Various'}}

@pytest.fixture
def expected_order_info():
    return {'Primates': ['Gibbon', 'Olive_baboon'],
            'Eulipotyphla': ['Hedgehog'],
            'Xenarthra': ['Armadillo'],
            'Cetacea': ['Blue_whale'],
            'Chiroptera': ['Megabat']}

@pytest.fixture
def order_reference():
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_order_info.json', 'r') as infile:
        order_reference = json.load(infile)
    return order_reference

def test_extract_species_names(sequences, expected_species_names):
    species_names = extract_species_names(sequences)

    assert species_names == expected_species_names, "Species_name_don't match"


def test_extract_taxonomic_order_info(taxa_reference, order_reference, sequences, expected_taxa_info, expected_order_info):
    species_names = extract_species_names(sequences)

    # Call the function
    taxanomic_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)

    # Assert the expected outcomes
    for species_name in species_names:
        assert taxanomic_info[species_name] == expected_taxa_info[species_name], "Taxonomic info does not match expected"
    
    for order, species in expected_order_info.items():
        assert order_info[order] == species, "Order info does not match expected"

def test_pseudorandom_sequence_selection_within_order(taxa_reference, order_reference, sequences):
    species_names = extract_species_names(sequences)

    traid_selected = pseudorandom_sequence_selection_within_order(species_names, taxa_reference, order_reference)
    
    ingroup1 = traid_selected['ingroup1']
    ingroup2 = traid_selected['ingroup2']
    outgroup = traid_selected['outgroup']

    taxanomic_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)
    
    assert taxanomic_info[ingroup1]['order'] == taxanomic_info[ingroup2]['order'], "Ingroup species order is not the same"
    assert taxanomic_info[outgroup]['order'] != taxanomic_info[ingroup1]['order'], "Out group and ingroup species has the same order"
    assert taxanomic_info[outgroup]['order'] != taxanomic_info[ingroup2]['order'], "Out group and ingroup species has the same order"

def test_pseudorandom_sequence_selection_between_order(taxa_reference, order_reference, sequences):
    species_names = extract_species_names(sequences)

    traid_selected = pseudorandom_sequence_selection_between_order(species_names, taxa_reference, order_reference)
    
    ingroup1 = traid_selected['ingroup1']
    ingroup2 = traid_selected['ingroup2']
    outgroup = traid_selected['outgroup']

    taxanomic_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)
    
    assert taxanomic_info[outgroup]['order'] != taxanomic_info[ingroup1]['order'], "Out group and ingroup species has the same order"
    assert taxanomic_info[outgroup]['order'] != taxanomic_info[ingroup2]['order'], "Out group and ingroup species has the same order"




def test_traids_selecting_binning():
    # Setup test data
    traid_data = {'traits_speies_name': {'ingroup1': 'Olive_baboon', 'ingroup2': 'Gibbon', 'outgroup': 'Blue_whale'},
        'traits_info': {
            'ingroup_jsd': 0.0005,
            'ens_difference': 0.0007,
            'shortest_ingroup_ens': 0.0003
        }
    }
    jsd_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 1)}
    ens_diff_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 1)}
    shortest_ens_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 1)}

    # Execute function
    updated_jsd_bins, updated_ens_diff_bins, updated_shortest_ens_bins, triad_added = traids_selecting_binning(traid_data, jsd_bins, ens_diff_bins, shortest_ens_bins)
    # Assert that triad was added to bins
    assert triad_added, "Triad should have been marked as added"
    assert updated_jsd_bins['bins'][0] is not None, "JSD bin should have data"
    assert updated_ens_diff_bins['bins'][0] is not None, "ENS difference bin should have data"
    assert updated_shortest_ens_bins['bins'][0] is not None, "Shortest ENS bin should have data"


def test_get_traits_info(sequences, taxa_reference, order_reference):
    species_names = extract_species_names(sequences)
    traits_species_names = pseudorandom_sequence_selection_within_order(species_names, taxa_reference, order_reference)

    traid_data = get_traits_info(traits_species_names, 'TN93', sequences)

    assert traid_data != None

