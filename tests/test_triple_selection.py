import pytest
from clock_project.genome_analysis.triple_selection import get_ingroup_pair, get_triples_info, get_relationship, get_outgroup_species, get_triples_info
from cogent3 import get_app

load_json_app = get_app("load_json")

@pytest.fixture
def sequence_path():
    return '/Users/gulugulu/Desktop/honours/data_local_2/sequences_filterd_550_threshold/ENSG00000265203.json'

@pytest.fixture
def result_lf():
    model_fit_result_path = '/Users/gulugulu/Desktop/honours/data_local_2/whole_gene_model_fitting_550_threshold/ENSG00000265203.json'
    load_json_app = get_app("load_json")
    result= load_json_app(model_fit_result_path)
    return result.lf

@pytest.fixture
def ens_tree(result_lf):
    ens_tree = result_lf.get_ens_tree()
    return ens_tree


@pytest.fixture
def sequences(sequence_path):
    sequences = load_json_app(sequence_path)
    return sequences

@pytest.fixture
def triples_species_name():
    return {'ingroup1': 'Arctic_ground_squirrel', 'ingroup2': 'Ryukyu_mouse', 'outgroup': 'Pika'}


@pytest.fixture
def expected_relationship():
    return {'sibling': ['Kangaroo_rat',
  'Ryukyu_mouse',
  'Rat',
  'Golden_Hamster',
  'Prairie_vole',
  'Northern_American_deer_mouse',
  'Lesser_Egyptian_jerboa',
  'Long_tailed_chinchilla',
  'Degu',
  'Guinea_Pig',
  'Arctic_ground_squirrel',
  'Squirrel',
  'Alpine_marmot',
  'Eurasian_red_squirrel'],
 'cousion': ['Rabbit', 'Pika']}

# @pytest.fixture
# def expected_Q_matrices():
#     return {'Gibbon': np.array([[-0.84964026,  0.57954512,  0.1837098 ,  0.08638534],
#         [ 1.16322517, -1.5475045 ,  0.15549139,  0.22878794],
#         [ 0.32121906,  0.15514015, -1.04933063,  0.57297143],
#         [ 0.12209357,  0.05652131,  0.71071505, -0.88932993]]),
#  'Rat': np.array([[-1.46189468,  0.97752474,  0.24307002,  0.24129992],
#         [ 0.72461046, -0.91704654,  0.11476656,  0.07766951],
#         [ 0.23491255,  0.39778288, -1.21807961,  0.58538418],
#         [ 0.02740988,  0.12992142,  0.49959346, -0.65692477]]),
#  'Sperm_whale': np.array([[-2.00330541,  1.49496865,  0.27115223,  0.23718453],
#         [ 1.7271533 , -2.13730183,  0.40176229,  0.00838624],
#         [ 0.169661  ,  0.27751602, -1.43616211,  0.98898508],
#         [ 0.6064174 ,  0.14670557,  0.60246588, -1.35558885]])}


@pytest.fixture
def triples_info_big_tree():
    return {'triples_info_big_tree': {'internal_node': 'Rodentia.3',
  'root': 'Glires',
  'jsd_dict': {'Ingroup_JSD': 0.0016731612268652274,
   'JSD_difference': {'Arctic_ground_squirrel': 0.002180861209883478,
    'Ryukyu_mouse': 0.004796972640505048}},
  'ens_difference': 0.16052010587199705,
  'shortest_ingroup_ens': 0.18970312896624109,
  'ens': {'Arctic_ground_squirrel': 0.18970312896624109,
   'Ryukyu_mouse': 0.35022323483823814,
   'Pika': 0.3893442432544122},
  'nuc_freqs_dict': {'Ingroup_nuc_freqs': {'Arctic_ground_squirrel': [0.11831442463533225,
     0.35170178282009723,
     0.3500810372771475,
     0.17990275526742303],
    'Ryukyu_mouse': [0.146677471636953,
     0.3557536466774716,
     0.31766612641815234,
     0.17990275526742303]},
   'internal_root_distribution': [0.14258068426038356,
    0.38155418452218853,
    0.11135943021404651,
    0.364505701003382],
   'root_distribution': [0.14437639788652204,
    0.37970043986338864,
    0.11241661954959037,
    0.3635065427004995]}}}

@pytest.fixture
def lca_parent_node():
    return 'Rodentia.3', 'Glires'


def test_get_ingroup_pair(sequences):
    ingroup_pair = get_ingroup_pair(sequences)
    all_speices = sequences.names
    assert set(ingroup_pair) <= set(all_speices), "ingroup pair is not species in the sequence."

def test_get_outgroup_species(triples_species_name, ens_tree, expected_relationship):
    ingroup_species_pair = [triples_species_name['ingroup1'], triples_species_name['ingroup2']]
    outgroup = get_outgroup_species(ingroup_species_pair, ens_tree)

    assert set(outgroup) <= set(expected_relationship['cousion']), "Outgroup species is not the sibling of lca."


def test_get_relationship(triples_species_name, expected_relationship, ens_tree):
    ingroup_species_pair = [triples_species_name['ingroup1'], triples_species_name['ingroup2']]
    relationship = get_relationship(ingroup_species_pair, ens_tree)

    assert relationship == expected_relationship, "Species_name_don't match"

def test_get_triples_info(triples_species_name, ens_tree, triples_info_big_tree, lca_parent_node, sequences, result_lf):
    lca_node, parent_node = lca_parent_node 
    traids_info = get_triples_info(triples_species_name, ens_tree, sequences, result_lf, lca_node, parent_node, iteration_count=2)
    assert traids_info[0]['triples_info_big_tree'] == triples_info_big_tree['triples_info_big_tree'], 'there is difference in triples_info.'



    

