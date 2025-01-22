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
    output_file = os.path.join(output_dir, f"{os.path.basename(file)}_counts.csv")
    pd.DataFrame(list(counts.items()), columns=["Sequence", "Count"]).to_csv(output_file, index=False)
    return output_file

def merge_intermediate_files(intermediate_files, output_file):
    """
    Merge all intermediate CSV files into a single consolidated DataFrame and save as CSV.
    """
    merged_counts = {}

    for file in intermediate_files:
        df = pd.read_csv(file)
        for _, row in df.iterrows():
            seq = row["Sequence"]
            count = row["Count"]
            merged_counts[seq] = merged_counts.get(seq, 0) + count

    # Convert merged counts to DataFrame
    merged_df = pd.DataFrame(list(merged_counts.items()), columns=["Sequence", "Count"])
    merged_df.sort_values(by="Count", ascending=False, inplace=True)
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
    file_list = ["sample1.fastq", "sample2.fastq"]  # List of input FASTQ files
    primer_fwd = "ATCG"  # Forward primer sequence
    primer_rev = "GCTA"  # Reverse primer sequence
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
