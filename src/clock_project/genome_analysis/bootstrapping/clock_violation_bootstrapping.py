from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_moltype, get_app
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from cogent3.app import evo
RATE_PARAM_UPPER = 50

def get_id(result):
    return result.source.unique_id

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


@define_app
def test_hypothesis_clock_model_N(aln: AlignedSeqsType, tree=None, opt_args=None) -> SerialisableType:
    outgroup_name = aln.info['triples_species_name']['outgroup']
    tree = make_tree(tip_names=aln.names)
    print(outgroup_name)
    outgroup_edge = [outgroup_name]

    model_kwargs = dict(
    tree=tree,
    opt_args=opt_args,
    # unique_trees=True,
    lf_args=dict(discrete_edges=[outgroup_edge]),
    optimise_motif_probs=True,
    )
    null = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            time_het = None,
            **model_kwargs,
        )
    alt = evo.model(
            "GN",
            name="GN-max-het",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            time_het = "max",
            **model_kwargs,
        )
    
    hyp = evo.hypothesis(null, alt, sequential=True)
    bootstrapper = evo.bootstrap(hyp, num_reps=100, parallel=True)
    result = bootstrapper(aln)    
    print('finish')
    return result

load_json_app = get_app("load_json")

def p_value(result):
    return sum(result.observed.LR <= null_lr for null_lr in result.null_dist) / len(result.null_dist)

clock_bootstrapper = test_hypothesis_clock_model_N()

