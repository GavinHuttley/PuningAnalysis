import pytest
import numpy as np
import json
from cogent3.util.deserialise import deserialise_object
from clock_project.genome_analysis.triad_selection_2 import get_ingroup_pair, get_triads_info, get_relationship, get_outgroup_species, get_Q_matrices, get_triads_info
from cogent3 import get_app
@pytest.fixture
def sequence_path():
    return '/Users/gulugulu/repos/PuningAnalysis/data/ensembl_ortholog_sequences/homologies_sequence_common_name/ENSG00000162620.json'

@pytest.fixture
def result_lf():
    model_fit_result_path = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result/ENSG00000162620.json'
    load_json_app = get_app("load_json")
    result_lf = load_json_app(model_fit_result_path)
    return result_lf

@pytest.fixture
def ens_tree(result_lf):
    ens_tree = result_lf.get_ens_tree()
    return ens_tree


@pytest.fixture
def sequences(sequence_path):
    sequences = deserialise_object(json.load(open(sequence_path, 'r')))
    return sequences

@pytest.fixture
def triads_species_name():
    return {'ingroup1': 'Vervet_AGM', 'ingroup2': 'Cat', 'outgroup': 'Armadillo'}

@pytest.fixture
def traids_species_name():
    return {'ingroup1': 'Vervet_AGM', 'ingroup2': 'Cat', 'outgroup': 'Armadillo'}

@pytest.fixture
def expected_relationship():
    return {'sibling': ['American_black_bear', 'Giant_panda', 'Ferret', 'Dog', 'Red_fox', 'Cat', 'Lion', 'Donkey', 'Vaquita',
  'Dolphin',
  'Blue_whale',
  'Goat',
  'Domestic_yak',
  'American_bison',
  'Siberian_musk_deer',
  'Yarkand_deer',
  'Pig',
  'Chacoan_peccary',
  'Arabian_camel',
  'Alpaca',
  'Greater_horseshoe_bat',
  'Kangaroo_rat',
  'Ryukyu_mouse',
  'Rat',
  'Golden_Hamster',
  'Prairie_vole',
  'Northern_American_deer_mouse',
  'Long_tailed_chinchilla',
  'Degu',
  'Guinea_Pig',
  'Squirrel',
  'Alpine_marmot',
  'Eurasian_red_squirrel',
  'Rabbit',
  'Chimpanzee',
  'Gorilla',
  'Sumatran_orangutan',
  'Gibbon',
  'Golden_snub_nosed_monkey',
  'Sooty_mangabey',
  'Drill',
  'Olive_baboon',
  'Macaque',
  'Vervet_AGM',
  'Panamanian_white_faced_capuchin',
  'Bolivian_squirrel_monkey',
  "Ma's_night_monkey",
  'White_tufted_ear_marmoset',
  'Tarsier',
  'Bushbaby',
  'Tree_Shrew'],
 'cousion': ['Armadillo']}

@pytest.fixture
def expected_Q_matrices():
    return {'Vervet_AGM': np.array([[-1.53172579e+00,  6.58813542e-01,  7.96540377e-01,
          7.63718671e-02],
        [ 1.49016212e+00, -1.46302728e+00, -2.61252689e-02,
         -1.00957960e-03],
        [ 5.48426316e-08,  5.36226704e-08, -1.29376043e+00,
          1.29376032e+00],
        [ 5.48833587e-08,  5.36424254e-08,  1.38790347e-04,
         -1.38898873e-04]]),
 'Cat': np.array([[-1.07690709e+00,  5.50676954e-01,  3.10509477e-01,
          2.15720659e-01],
        [ 8.82658154e-01, -1.74376389e+00,  8.75257325e-01,
         -1.41515921e-02],
        [-4.47929147e-03,  2.56654788e-01, -1.21975245e+00,
          9.67576956e-01],
        [ 3.45084021e-05,  4.54342972e-05,  6.05028939e-01,
         -6.05108881e-01]]),
 'Armadillo': np.array([[-1.80901024e+00,  1.19022820e+00,  4.94123117e-01,
          1.24658924e-01],
        [ 2.32262640e+00, -2.52365530e+00,  1.99985438e-01,
          1.04346846e-03],
        [ 1.59853799e-01, -6.33027107e-03, -9.65703113e-01,
          8.12179586e-01],
        [ 2.52027188e-01,  7.03450723e-03,  1.99922168e+00,
         -2.25828337e+00]])}




@pytest.fixture
def expected_triads_info():
    return {'triads_species_names': {'ingroup1': 'Vervet_AGM',
  'ingroup2': 'Cat',
  'outgroup': 'Armadillo'},
 'triads_info_big_tree': {'internal_node': 'Boreoeutheria',
  'root': 'Eutheria.2',
  'ingroup_jsd': 0.0035286444750179946,
  'ens_difference': 0.034038717242672614,
  'shortest_ingroup_ens': 0.04412841189885737,
  'ens': {'Vervet_AGM': 0.04412841189885737,
   'Cat': 0.07816712914152998,
   'Armadillo': 0.1420254581405398},
  'internal_root_distribution': [0.4039741240439267,
   0.17448172346622864,
   0.3108684352076752,
   0.11067571728216971],
  'root_distribution': [0.4155761298210775,
   0.17448172329851894,
   0.299493919991724,
   0.11044822688867967],
  'nabla_values': {'Vervet_AGM': 0.5501733102511572,
   'Cat': 0.39590857339503965,
   'Armadillo': 0.25066639863458745}}}

@pytest.fixture
def lca_parent_node():
    return 'Boreoeutheria', 'Eutheria.2'


def test_get_ingroup_pair(sequences):
    ingroup_pair = get_ingroup_pair(sequences)
    all_speices = sequences.names
    assert set(ingroup_pair) <= set(all_speices), "ingroup pair is not species in the sequence."

def test_get_relationship(triads_species_name, expected_relationship, ens_tree):
    ingroup_species_pair = [triads_species_name['ingroup1'], triads_species_name['ingroup2']]
    relationship = get_relationship(ingroup_species_pair, ens_tree)

    assert relationship == expected_relationship, "Species_name_don't match"


def test_get_outgroup_species(triads_species_name, ens_tree, expected_relationship):
    ingroup_species_pair = [triads_species_name['ingroup1'], triads_species_name['ingroup2']]
    outgroup = get_outgroup_species(ingroup_species_pair, ens_tree)

    assert set(outgroup) <= set(expected_relationship['cousion']), "Outgroup species is not the sibling of lca."

def test_get_Q_matrix(triads_species_name, ens_tree, expected_Q_matrices, lca_parent_node, result_lf):
    lca_node, parent_node = lca_parent_node 
    sub_tree = ens_tree.get_sub_tree(triads_species_name.values())
    ens_dict = {n: sub_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triads_species_name.values()}
    Q_matrices = get_Q_matrices(triads_species_name, lca_node, parent_node, ens_tree, result_lf, ens_dict)
    assert Q_matrices.all() == expected_Q_matrices.all(), 'rate matrix is not as expected'


def test_get_triads_info(triads_species_name, ens_tree, expected_triads_info, lca_parent_node, sequences, result_lf):
    lca_node, parent_node = lca_parent_node 
    traids_info = get_triads_info(triads_species_name, ens_tree, sequences, result_lf, lca_node, parent_node)
    assert traids_info[0] == expected_triads_info, 'there is difference in triads_info.'



    

