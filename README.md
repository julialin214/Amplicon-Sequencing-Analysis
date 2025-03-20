# Amplicon Sequencing Analysis 

### `sequence_counter.py`
Counts unique sequences across multiple FASTQ files. The script processes a list of FASTQ inputs, extracts sequences based on specified 5' and 3' constant regions, and filters them by a defined variable region length range. It outputs a CSV file with the counts of each unique sequence across all input files.

### `RBS_base_comp.py`
Analyzes the base composition of the ribosome binding site (RBS) in DNA sequences extracted from a CSV file. It extracts and quantifies nucleotide composition in a 6-nucleotide region 6 nucleotides from the 3' end (assumed to be the start codon) 

Input File Format
- The CSV file should have sample names in the first row.
- Each column represents a sample, with DNA sequences listed below the sample name.
- Example Input (CSV format):

Sample1	Sample2	Sample3

ACGTG...	TTAGC...	GGTCA...

TATAG...	ATCGG...	GCTAT...
