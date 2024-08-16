from cogent3 import get_app, open_data_store
from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_moltype
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from cogent3.app import evo

@register_model("nucleotide")
def GSN(**kwargs):
    """A General Stationary Nucleotide substitution model instance."""
    kwargs["optimise_motif_probs"] = kwargs.get("optimise_motif_probs", True)
    kwargs["name"] = kwargs.get("name", "GSN")
    return GeneralStationary(get_moltype("dna").alphabet, **kwargs)

def get_param_rules_upper_limit(model_name, upper):
    """rules to set the upper value for rate matrix terms"""
    from cogent3 import get_model

    sm = get_model(model_name)
    return [{"par_name": par_name, "upper": upper} for par_name in sm.get_param_list()]


RATE_PARAM_UPPER = 50

@define_app
def test_hypothesis_model_bootstrapper(aln: AlignedSeqsType, tree=None, opt_args=None) -> SerialisableType:
    outgroup_name = aln.info['triads_species_name']['outgroup']
    tree = make_tree(tip_names=aln.names)
    
    outgroup_edge = [outgroup_name]

    model_kwargs = dict(
    tree=tree,
    opt_args=opt_args,
    # unique_trees=True,
    lf_args=dict(discrete_edges=[outgroup_edge]),
    optimise_motif_probs=True,
    )
    null = evo.model(
            "GSN",
            param_rules=get_param_rules_upper_limit("GSN", RATE_PARAM_UPPER),
            **model_kwargs,
        )
    alt = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            **model_kwargs,
        )
    
    hyp = evo.hypothesis(null, alt, sequential=True)
    bootstrapper = evo.bootstrap(hyp, num_reps=1, parallel=True)
    result = bootstrapper(aln)    
    return result

def p_value(result):
    return sum(result.observed.LR <= null_lr for null_lr in result.null_dist) / len(result.null_dist)



bootstrapper = test_hypothesis_model_bootstrapper()


# aln_dir_new = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset_info_added'

# input_data_store = open_data_store(aln_dir_new, suffix= 'json')
# bootstrap_results = list(path_bootstrapper.as_completed(input_data_store))









