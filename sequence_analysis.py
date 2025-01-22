from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool
from functools import partial
import os

def process_file(file, primer_fwd, primer_rev, min_len, max_len):
    """
    Process a single FASTQ file to extract and count regions of interest.
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
    
    return file, counts


def process_files_parallel(file_list, primer_fwd, primer_rev, min_len, max_len, num_processes):
    """
    Process multiple files in parallel using multiprocessing.
    """
    # Partial function to pass constant parameters
    partial_process = partial(process_file, primer_fwd=primer_fwd, primer_rev=primer_rev, min_len=min_len, max_len=max_len)
    
    with Pool(num_processes) as pool:
        results = pool.map(partial_process, file_list)
    
    return results

def merge_results(results):
    """
    Merge counts from all files into a single DataFrame.
    """
    merged_data = {}
    
    for file, counts in results:
        sample_name = os.path.basename(file).split(".")[0]
        
        for region, count in counts.items():
            if region not in merged_data:
                merged_data[region] = {}
            merged_data[region][sample_name] = count
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(merged_data, orient="index").fillna(0)
    df.index.name = "Sequence"
    
    return df

def main():
    """
    Main function to execute the program.
    """
    # User input
    #file_list = ["sample1.fastq", "sample2.fastq"]  # List of input FASTQ files
    file_list = ["sample1.fastq"]  # List of input FASTQ files
    primer_fwd = "ATCG"  # Forward primer sequence
    primer_rev = "GCTA"  # Reverse primer sequence
    min_len = 0         # Minimum length of region of interest
    max_len = 100        # Maximum length of region of interest
    output_file = "sequence_counts.csv"
    num_processes = 4    # Number of parallel processes
    
    # Process files in parallel
    print("Processing files...")
    results = process_files_parallel(file_list, primer_fwd, primer_rev, min_len, max_len, num_processes)
    
    # Merge results
    print("Merging results...")
    df = merge_results(results)
    
    # Export to CSV
    print(f"Exporting results to {output_file}...")
    df.to_csv(output_file)
    print("Done!")

if __name__ == "__main__":
    main()
