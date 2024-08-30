import pytest
from cogent3 import get_app, open_data_store
import os
import json

load_json_app = get_app("load_json")

@pytest.fixture
def rep_subset_dir():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_aln_representative_subset'
@pytest.fixture
def resource_dir():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_alignments_550_threshold'

@pytest.fixture
def info_dir():
    return '/Users/gulugulu/Desktop/honours/data_local_2/triples_550_threshold'


def test_correct_file_added(rep_subset_dir, resource_dir):
    result_input_sdtore = open_data_store(rep_subset_dir, mode="r", suffix="json")
    for result in result_input_sdtore:
        gene_name_identifier = result.unique_id.split('.')[0]
        gene_name, identifier = gene_name_identifier.split('_')
        res = load_json_app(result)
        resource_path = os.path.join(resource_dir, gene_name, f"{identifier}.json")
        resource = load_json_app(resource_path)
        assert res == resource

def test_species_name_exist(rep_subset_dir, info_dir):
    result_input_sdtore = open_data_store(rep_subset_dir, mode="r", suffix="json")
    for result in result_input_sdtore:
        gene_name_identifier = result.unique_id.split('.')[0]
        gene_name, identifier = gene_name_identifier.split('_')
        res = load_json_app(result)
        info_path = os.path.join(info_dir, gene_name, "triples_species_names_dict.json")
        info_re = json.load(open(info_path, 'r'))
        triples_species_name = info_re[identifier]
        assert res.info['triples_species_name'] == triples_species_name

        

