import os
import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial

def process_file(file, primer_fwd, primer_rev, min_len, max_len, output_dir):
    """
    Process a single FASTQ file to extract and count regions of interest.
    Write intermediate results to a CSV file.
    """
    counts = {}

    for record in SeqIO.parse(file, "fastq"):
        seq = str(record.seq)

        # Locate primers
        start = seq.find(primer_fwd)
        if start == -1:  # Forward primer not found
            continue
        start += len(primer_fwd)  # Adjust start position to the end of the forward primer

        end = seq.find(primer_rev, start)
        if end == -1:  # Reverse primer not found
            continue

        # Extract sequence between primers
        region = seq[start:end]

        # Check length constraints
        if min_len <= len(region) <= max_len:
            counts[region] = counts.get(region, 0) + 1

    # Write intermediate results to a CSV file
    sample_name = os.path.basename(file).split(".")[0]
    output_file = os.path.join(output_dir, f"{sample_name}_counts.csv")
    pd.DataFrame(list(counts.items()), columns=["Sequence", sample_name]).to_csv(output_file, index=False)
    return output_file

def merge_intermediate_files(intermediate_files, output_file):
    """
    Merge all intermediate CSV files into a single DataFrame with separate columns for each file.
    """
    merged_df = None

    for file in intermediate_files:
        # Read intermediate file
        df = pd.read_csv(file)
        sample_name = os.path.basename(file).split("_counts.csv")[0]

        if merged_df is None:
            # Initialize the merged DataFrame with the first file
            merged_df = df
        else:
            # Merge on the "Sequence" column, keeping separate columns for counts
            merged_df = pd.merge(merged_df, df, on="Sequence", how="outer")

    # Fill NaN with 0 for sequences not present in all files
    merged_df.fillna(0, inplace=True)

    # Save the merged DataFrame
    merged_df.to_csv(output_file, index=False)

def process_files_parallel(file_list, primer_fwd, primer_rev, min_len, max_len, num_processes, output_dir):
    """
    Process multiple files in parallel using multiprocessing and save intermediate results to disk.
    """
    os.makedirs(output_dir, exist_ok=True)
    partial_process = partial(process_file, primer_fwd=primer_fwd, primer_rev=primer_rev, min_len=min_len, max_len=max_len, output_dir=output_dir)

    with Pool(num_processes) as pool:
        intermediate_files = pool.map(partial_process, file_list)

    return intermediate_files

def main():
    """
    Main function to execute the program.
    """
    # User input
    file_list = ["JLSX2b.dedup.merge.fastq"]  # List of input FASTQ files
    primer_fwd = "CGTAGTCGTAGCTGATCGAC"  # Forward primer sequence
    primer_rev = "ATGTCTCTAAGTACTGAA"  # Reverse primer sequence
    min_len = 0          # Minimum length of region of interest
    max_len = 100        # Maximum length of region of interest
    output_dir = "intermediate_counts"  # Directory for intermediate CSV files
    final_output = "merged_sequence_counts.csv"  # Final merged output file
    num_processes = 4     # Number of parallel processes

    # Process files in parallel
    print("Processing files in parallel...")
    intermediate_files = process_files_parallel(file_list, primer_fwd, primer_rev, min_len, max_len, num_processes, output_dir)

    # Merge intermediate results
    print("Merging intermediate results...")
    merge_intermediate_files(intermediate_files, final_output)

    print(f"Final results saved to {final_output}")

if __name__ == "__main__":
    main()
