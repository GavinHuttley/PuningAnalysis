import pytest
from pathlib import Path
from cogent3 import open_data_store, get_app

# Adjust this to point to your test bootstrap result file
toc_result_path = '~/clock/mammal_orthologs_hsap_1/toc_result'

loader_json = get_app('load_json')
toc_result_dstore = open_data_store(toc_result_path, suffix = 'json')
tabulate = get_app("tabulate_stats")


@pytest.fixture(scope="module")
def bootstrap_result():
    """Load the bootstrap result"""
    dstore = open_data_store(toc_result_path, suffix="json")
    results = list(dstore)
    assert len(results) > 0, "No bootstrap result files found"
    return loader_json(results[0])

def test_hypothesis_model_names(bootstrap_result):
    """Test if the null and alt model names are correct"""
    hyp_result = bootstrap_result

    assert hasattr(hyp_result, "observed"), "Result missing observed attribute"

    null_model = hyp_result.observed.null
    alt_model = hyp_result.observed.alt

    assert null_model.name == "clock", f"Expected null model name 'clock', got {null_model.name}"
    assert alt_model.name == "no-clock", f"Expected alt model name 'no-clock', got {alt_model.name}"

def test_bootstrap_replicates(bootstrap_result):
    """Test if the correct number of bootstrap replicates were run"""
    bootstraps = bootstrap_result.null_dist
    num_reps = len(bootstraps)
    assert num_reps == 100, "Bootstrap replicates incorrect"


def test_bootstrap_results_valid(bootstrap_result):
    """Test if bootstrap log-likelihood values look reasonable"""
    bootstraps_keys = bootstrap_result.keys()
    for i in bootstraps_keys:
        assert isinstance(bootstrap_result[i].LR, float), "LR value is not float"
        assert isinstance(bootstrap_result[i].null.lnL, float), "lnL value is not float"
        assert isinstance(bootstrap_result[i].alt.lnL, float), "lnL value is not float"

def test_bootstrap_nfp_diff(bootstrap_result):
    """Test if the correct number of bootstrap replicates were run"""
    df = bootstrap_result.observed.df
    assert df == 1, f"Expected 1 degree of freedom, got {df}"

def test_model_args(bootstrap_result):
    """Test if the model arguments are correctly set"""
    null_model = bootstrap_result.observed.null
    alt_model = bootstrap_result.observed.alt

    null_stats = tabulate(null_model)
    alt_stats = tabulate(alt_model)

    #expect mix discrete continous time model to have 'global params' and 'motif motif2 params'
    assert 'global params' in null_stats, "Null model args missing 'global params'"
    assert 'global params' in alt_stats, "Alt model args missing 'global params'"
    assert 'motif motif2 params' in alt_stats, "Alt model args missing 'motif params'"
    assert 'motif motif2 params' in null_stats, "Null model args missing 'motif params'"
    #for the clock, no 'edge params' in the null model, only one rate matrix/edge length
    assert 'edge params' in alt_stats, "Alt model args missing 'edge params'"
    assert 'edge params' not in null_stats, "Null model should not have 'edge params'"







    




    



# def test_serialization_roundtrip(bootstrap_result, tmp_path):
#     """Test that bootstrap result can be re-serialized/deserialized correctly"""
#     temp_file = tmp_path / "temp_result.json"
#     bootstrap_result.write(temp_file)

#     # re-read and deserialize
#     new_result = deserialise_object(temp_file.read_text())
#     assert new_result.model.null.name == "clock"
#     assert new_result.model.alt.name == "no-clock"
#     assert len(new_result.bootstraps) == len(bootstrap_result.bootstraps)
