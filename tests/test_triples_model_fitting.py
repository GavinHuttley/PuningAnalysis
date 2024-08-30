import pytest
from clock_project.genome_analysis.triples_model_fitting import check_process, triples_model_fitting
from cogent3 import get_app, open_data_store

@pytest.fixture
def aln_path():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_alignments_550_threshold/ENSG00000009307'

@pytest.fixture
def aln_path():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_alignments_550_threshold/ENSG00000009307'


@pytest.fixture
def triples_species_names():
    return {'ingroup1': 'Elephant', 'ingroup2': 'Prairie_vole', 'outgroup': 'Shrew'}

@pytest.fixture
def result_lf_dir():
    return '/Users/gulugulu/Desktop/honours/data_local_2/whole_gene_model_fitting_550_threshold'

@pytest.fixture
def triples_info_dir():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_550_threshold'

def test_check_process(aln_path, result_lf_dir, triples_info_dir):
    triples_speccies_infos, result, triples_aln = check_process(aln_path, result_lf_dir, triples_info_dir)
    triples_alignment_paths = open_data_store(aln_path, mode="r", suffix="json")
    identifier_number = triples_alignment_paths[0].unique_id.split('.')[0]

    assert set(triples_speccies_infos[identifier_number].values()) == set(triples_aln.names)
    assert result.source.split('.')[0] == 'ENSG00000009307'

def test_triples_model_fitting(aln_path, result_lf_dir, triples_info_dir, triples_species_names):
    triples_species_infos, result, triples_aln = check_process(aln_path, result_lf_dir, triples_info_dir)
    triples_alignment_paths = open_data_store(aln_path, mode="r", suffix="json")
    identifier_number = triples_alignment_paths[0].unique_id.split('.')[0]
    ens_tree = result.lf.get_ens_tree()
    model_fitting_result = triples_model_fitting(triples_aln, triples_species_names, ens_tree)

    assert model_fitting_result.alignment == triples_aln
    assert triples_species_infos[identifier_number] == triples_species_names
    assert set(model_fitting_result.tree.get_tip_names()) == set(triples_aln.names)
    assert model_fitting_result.tree.get_newick() == ens_tree.get_sub_tree(triples_species_names.values()).get_newick()