import json
from yapeng_check_BV import get_bounds_violation
from cogent3.util.deserialise import deserialise_object
import numpy as np
import click
import os
import glob

bounary_violation_function = get_bounds_violation()


base_dir = '/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_model_fitting_info/ENSG00000116285'

model_fitting_result_dir = os.path.join(base_dir, 'model_fitting_result')

model_fitting_results_path = glob.glob(os.path.join(model_fitting_result_dir, '*.json'))

bounary_violation_check = []

for path in model_fitting_result_dir:
    identifier_number = os.path.basename(path).rsplit('.', 1)[0]
    reuslt_serialised = json.load(open(path, 'r'))
    model_result = deserialise_object()
    bounary_violarion = bounary_violation_function.main(model_result)



