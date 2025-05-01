import dataclasses

from cogent3.app.composable import define_app
from cogent3.app.result import model_result



@dataclasses.dataclass(slots=True)
class ParamRules:
    source: str
    params: list[dict]


@dataclasses.dataclass(slots=True)
class BoundsViolation:
    source: str
    vio: list[dict]




def load_param_values(model_result:model_result) -> ParamRules:
    """get non-topology param values"""
    return ParamRules(
        source=model_result.source, params=model_result.lf.get_param_rules()
    )


class get_bounds_violation:
    """Check if there are any rate params proximity to the bounds as 1e-10.
    This value is important as two clusters of fits are split by the value."""

    exclude_params = ("length", "mprobs")  # Tuple should be enclosed in parentheses

    def __init__(self) -> None:
        pass

    def main(self, model_result) -> BoundsViolation:
        params = load_param_values(model_result)
        vio = []
        list_of_params = params.params
        for param in list_of_params:
            if param["par_name"] not in self.exclude_params:
                if (abs(param["init"] - param["lower"]) <= 1e-12) or (
                    abs(param["init"] - param["upper"]) <= 1e-12):
                    vio.append(param)
        return BoundsViolation(source=params.source, vio=vio)


