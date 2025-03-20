# Analyze Base Composition of Ribosome Binding Site 
#
# Installation:
# pip install pandas
#
# Parameters:
# - input_csv: The name of the input CSV file containing sequencing data.
# - output_dir: The directory where output files will be saved.
#
# Running the script:
# python RBS_base_comp.py

import pandas as pd
import os

def count_base_comp(input_csv):
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Get the sample names from the first row
    sample_names = df.columns.tolist()
    
    # Create output directory if it doesn't exist
    output_dir = "RBS_output_files"
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each sample column
    for sample in sample_names:
        sequences = df[sample].dropna().tolist()  # Remove NaN values
        
        # Prepare output data
        output_data = []
        
        for seq in sequences[1:]:  # Skip the first row (sample name)
            if len(seq) < 12:
                continue  # Skip sequences too short for analysis
            
            last_six = seq[-12:-6]  # Extract the target 6 nucleotides
            a_count = last_six.count('A')
            g_count = last_six.count('G')
            c_count = last_six.count('C')
            t_count = last_six.count('T')
            
            purine_count = a_count + g_count
            purine_fraction = purine_count / 6  # Fraction of purines in the 6-nt region
            
            output_data.append([seq, a_count, g_count, c_count, t_count, purine_count, purine_fraction])
        
        # Convert to DataFrame
        output_df = pd.DataFrame(output_data, columns=[sample, 'A_count', 'G_count', 'C_count', 'T_count', 'Purine_count', 'Purine_fraction'])
        
        # Save to CSV file
        output_file = os.path.join(output_dir, f"{sample}_RBS_base_comp.csv")
        output_df.to_csv(output_file, index=False)
        print(f"Saved output for {sample} to {output_file}")

# Run the analysis
input_csv = "Tt 2-1 R1 ratio > 50.csv"
count_base_comp(input_csv)
