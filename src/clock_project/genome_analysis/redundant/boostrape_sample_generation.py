import os
import shutil
import random

def generate_bootstrap_samples(source_dir, output_dir, num_samples=100):
    # Create the main output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all file names in the source directory
    file_names = [os.path.join(source_dir, f) for f in os.listdir(source_dir) if f.endswith('.json')]
    num_files = len(file_names)
    
    # Generate each bootstrap sample
    for i in range(num_samples):
        sample_dir = os.path.join(output_dir, f"bootstrap_sample_{i+1}")
        os.makedirs(sample_dir, exist_ok=True)
        
        # Sample with replacement
        sampled_files = random.choices(file_names, k=num_files)
        
        # Copy files to the new directory
        file_count = {}
        for file_path in sampled_files:
            base_name = os.path.basename(file_path)
            if base_name in file_count:
                file_count[base_name] += 1
                new_name = f"{os.path.splitext(base_name)[0]}_{file_count[base_name]}{os.path.splitext(base_name)[1]}"
            else:
                file_count[base_name] = 1
                new_name = base_name
            
            shutil.copy(file_path, os.path.join(sample_dir, new_name))

# Set your directories here
source_dir = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset'
output_dir = '/Users/gulugulu/Desktop/honours/data_local/bootstrap_samples'

# Call the function to generate bootstrap samples
generate_bootstrap_samples(source_dir, output_dir)
