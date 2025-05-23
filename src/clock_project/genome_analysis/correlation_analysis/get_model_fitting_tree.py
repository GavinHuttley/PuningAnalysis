from ensembl_tui._species import Species
from cogent3 import load_tree
import click
import multiprocessing
import json
from cogent3 import get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3.app.typing import SeqsCollectionType, SerialisableType



with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/genome_information/common_names_mapping.json',  'r') as common_name_infile:
    common_names_mapping = json.load(common_name_infile)

tree_path = '/Users/gulugulu/repos/PuningAnalysis/data/dataset2_ensemble_trees/raw_data/tree_nh_file/vertebrates_species-tree_Ensembl.nh'

def decorate_vertebrate_tree(tree_path, common_names_mapping):
    """
    Decorate the vertebrate tree by reassigning tip names based on predefined mappings.

    Args:
        tree_path (str): Path to the Newick formatted tree file.
        common_names_mapping (dict): A dictionary for additional name mapping.

    Returns:
        Tree: The decorated tree object with updated tip names.
    """
    # Load the vertebrate tree from a file
    tree_vertebrates = load_tree(tree_path, format=None, underscore_unmunge=False)

    # Generate database prefixes from species names
    def make_db_prefixes():
        return [n.lower().replace(" ", "_") for n in Species.get_species_names()]

    db_prefixes = make_db_prefixes()

    # Find the database prefix for a given species name
    def find_db_prefix(name):
        for db_prefix in db_prefixes:
            if name.startswith(db_prefix):
                return db_prefix
        return None

    db_name_dict = {}
    probs = []

    # Create a mapping of tree tip names to database prefixes
    for tip_name in tree_vertebrates.get_tip_names():
        db_name = find_db_prefix(tip_name.lower())
        if db_name is None:
            db_name_dict[tip_name] = tip_name
            probs.append(tip_name)
        else:
            db_name_dict[tip_name] = db_name

    # Special handling for Mus musculus strains
    for key in db_name_dict.keys():
        if key.split('_')[0:2] == ['Mus', 'musculus']:
            if key != 'Mus_musculus_reference_CL57BL6_strain':
                db_name_dict[key] = 'mus_musculus_strain'
            else:
                db_name_dict[key] = 'mus_musculus'

    # Special handling for Sus scrofa strains
    for key in db_name_dict.keys():
        if key.split('_')[0:2] == ['Sus', 'scrofa']:
            if key != 'Sus_scrofa_reference_breed':
                db_name_dict[key] = 'sus_scrofa_strain'
            else:
                db_name_dict[key] = 'sus_scrofa'

    # Reassign names in the tree using the db_name_dict and common_names_mapping
    tree_vertebrates.reassign_names(db_name_dict)
    tree_vertebrates.reassign_names(common_names_mapping)

    return tree_vertebrates

decorated_vertebrates_tree = decorate_vertebrate_tree(tree_path, common_names_mapping)

loader = get_app("load_aligned", format="fasta", moltype="dna")

@define_app
def model_fitting_tree(alignment: SeqsCollectionType) -> SerialisableType:
    outgroup_species = ['aquila_chrysaetos_chrysaetos']
    species_names = alignment.names
    tip_names = species_names + outgroup_species
    tree_sub = decorated_vertebrates_tree.get_sub_tree(tip_names)
    tree_topology = tree_sub.get_sub_tree(species_names)
    model = get_app("model", "GN", tree=tree_topology, time_het='max', optimise_motif_probs=True)
    print('start')
    model_fitting_result = model(alignment)
    print('end')
    return model_fitting_result

@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('-o', '--output_dir', default='/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result', help='Output directory for JSON files')
@click.option('-n', '--num_processes', default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process')
def main(input, num_processes, output_dir, limit):
    """
    Fit model for sequences in different paths.
    
    INPUT: Path to the directory containing sequence files or a JSON file listing the paths.
    """

    input_aln_store = open_data_store(input, suffix= 'fa', mode="r", limit=limit)
    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    model_fitting_app = model_fitting_tree()
    writer_json = get_app("write_json", data_store=out_dstore)
    process = loader + model_fitting_app + writer_json

    process.apply_to(input_aln_store.completed, parallel=True, par_kw=dict(max_workers=num_processes))

if __name__ == "__main__":
    main()