import pandas as pd
import os
import argparse
import gzip
from Bio.Seq import Seq
from Bio import motifs
from collections import Counter

################################################################
# CLI 
arg_description = '''Script to identify random U6 barcodes 
longer than 20 bps'''

parser = argparse.ArgumentParser(description=arg_description)
parser.add_argument('input_dir', type=str, help='''Path to 
                    directory containing input fastqs''')

args = parser.parse_args()

input_dir = args.input_dir

################################################################
# Get list of Read 1 and Read 2 fastq files
print('Making and sorting lists of fastqs...')
files = os.listdir(input_dir)

R1_fastqs = [fastq for fastq in files if fastq.startswith('UBL')
             and fastq.endswith('R1_001.fastq.gz')]
R2_fastqs = [fastq for fastq in files if fastq.startswith('UBL')
             and fastq.endswith('R2_001.fastq.gz')]
# Sort lists so that file names are sorted by UBL numbers
R1_fastqs = sorted(R1_fastqs, key=lambda x: int(x.split('_')[0][3:]))
R2_fastqs = sorted(R2_fastqs, key=lambda x: int(x.split('_')[0][3:]))

################################################################
## Define the sequences flanking the barcodes
# Read 1
R1_left_flank = "TTGCTAGGACCGGCCTTAAAGCCGTACG"
R1_right_flank = "ACTAGTCTGGCCGTCGTTTTACAGTGAC"
# Read 2
R2_left_flank = "GTCACTGTAAAACGACGGCCAGACTAGT"
R2_right_flank = "CGTACGGCTTTAAGGCCGGTCCTAGCAA"

################################################################
## Identify barcodes, calculate average barcode length for each
## library and create an output dataframe of positive clones
## that have barcode > 45 bp. Expected random barcodes are
## around 50 bp. 
############################
# Define functions needed
# Sequence pattern matching function that allows mismatches
def find_with_mismatches(sequence, pattern, max_mismatches=2):
    for i in range(len(sequence) - len(pattern) + 1):
        mismatches = sum(1 for a, b in zip(sequence[i:i+len(pattern)], pattern) if a != b)
        if mismatches <= max_mismatches:
            return i
    return -1
###########################
# Function to extract barcode sequence from matched patterns
def extract_barcode(sequence, left_flank, right_flank):
    left_pos = find_with_mismatches(sequence, left_flank)
    right_pos = find_with_mismatches(sequence, right_flank)
    if left_pos != -1 and right_pos != -1:
        barcode = sequence[left_pos + len(left_flank):right_pos]
        return barcode
    else:
        return ""
###########################
# Function to calculate average barcode length
def calculate_avg_barcode_length(barcode_list):
    total_length = sum(len(sequence) for sequence in barcode_list)
    if barcode_list:
        average_length = total_length / len(barcode_list)
        return average_length
    else:
        return 0
###########################
# Function for creating barcode lists
def create_barcode_lists(R1_fastq, R2_fastq, R1_left_flank=R1_left_flank, 
                         R1_right_flank=R1_right_flank, 
                         R2_left_flank=R2_left_flank, 
                         R2_right_flank=R2_right_flank, input_dir=input_dir):
    with gzip.open(input_dir+'/'+R1_fastq, 'rt') as R1, gzip.open(input_dir+'/'+R2_fastq, 'rt') as R2:
        # Read first line (Header)
        R1_1 = R1.readline()
        R2_1 = R2.readline()

        # Create empty lists to store barcodes
        R1_barcode_list = []
        R2_barcode_list = []

        while (R1_1):
            # Read line 2 (sequence)
            R1_2 = R1.readline()
            R2_2 = R2.readline()
            # Read line 3 (+ seperator)
            R1_3 = R1.readline()
            R2_3 = R2.readline()
            # Read line 4 (quality score)
            R1_4 = R1.readline()
            R2_4 = R2.readline()

            # Extract barcode sequence
            R1_barcode = extract_barcode(R1_2, R1_left_flank, R1_right_flank)
            R1_barcode_list.append(R1_barcode)
            R2_barcode = extract_barcode(R2_2, R2_left_flank, R2_right_flank)
            R2_barcode_list.append(R2_barcode)

            # Read line 1 (header) for next cycle
            R1_1 = R1.readline()
            R2_1 = R2.readline()

        # Close File
        R1.close()
        R2.close()
        # Filter lists for empty barcodes
        R1_barcode_list = list(filter(None, R1_barcode_list))
        R2_barcode_list = list(filter(None, R2_barcode_list))
        # Return
        return R1_barcode_list, R2_barcode_list
# Create empty final results dataframe
positive_libs = pd.DataFrame(columns = ['R1_library_ID', 
                                     'R1_avg_barcode_length', 
                                     'R2_library_ID', 
                                     'R2_avg_barcode_length'])
###########################
## Loop to get results
print('Extracting barcode sequences and calculating average barcode lengths...')
# Read in files
for R1_fastq, R2_fastq in zip(R1_fastqs, R2_fastqs):
        R1_barcode_list, R2_barcode_list = create_barcode_lists(R1_fastq, R2_fastq)

        # Calculate average barcode lengths
        R1_avg_barcode_length = calculate_avg_barcode_length(R1_barcode_list)
        R2_avg_barcode_length = calculate_avg_barcode_length(R2_barcode_list)

        # Get library IDs
        R1_id = R1_fastq.split('_')[0]
        R2_id = R2_fastq.split('_')[0]

        # Make dataframe for results
        results = pd.DataFrame([{'R1_library_ID':R1_id, 
                                 'R1_avg_barcode_length':R1_avg_barcode_length,
                                 'R2_library_ID':R2_id,
                                 'R2_avg_barcode_length':R2_avg_barcode_length}])
        # Concatenate results with final results dataframe
        positive_libs = pd.concat([positive_libs, results])

###########################################
# Filter final results for barcodes > 45 bp
positive_libs = positive_libs[positive_libs['R1_avg_barcode_length']>=45]
positive_libs = positive_libs[positive_libs['R2_avg_barcode_length']>=45]
# Save final results as csv
outfile_name = 'positive_UBL_barcode_libs.csv'
outfile_save = input_dir+'/'+outfile_name
positive_libs.to_csv(outfile_save)

print(f'Results are saved here:{outfile_save}')

################################################################
## Get consensus barcode sequences for positive libraries.
## And create weblogo images for each consensus.
print('Determining Consensus sequences...')
#####################
## Function to only keep barcode sequences that have the most
## commonly occuring sequence length.
def common_sequence_length_filter(barcode_sequences):
    sequence_lengths = [len(seq) for seq in barcode_sequences]
    length_counts = Counter(sequence_lengths)
    most_common_length = length_counts.most_common(1)[0][0]
    filtered_sequences = [seq for seq in barcode_sequences if 
                          len(seq) == most_common_length]
    return filtered_sequences
#####################
# Create directory to store weblogo images
weblogo_path = f'{input_dir}/weblogo_consensus'
os.makedirs(weblogo_path, exist_ok=True)
####################
# Make new list of fastqs of positive libraries
R1_positive_ids = positive_libs['R1_library_ID'].to_list()
R2_positive_ids = positive_libs['R2_library_ID'].to_list()

R1_fastqs = []
for fastq in R1_fastqs:
    UBL_number = file.split('_')[0]
    if UBL_number in R1_positive_ids:
        R1_fastqs.append(fastq)

R2_fastqs = []
for fastq in R2_fastqs:
    UBL_number = file.split('_')[0]
    if UBL_number in R2_positive_ids:
        R2_fastqs.append(fastq)
######################
# Create results data frame
consensus_df = pd.DataFrame(columns = ['R1_library_ID', 
                                     'R1_consensus', 
                                     'R2_library_ID', 
                                     'R2_consensus'])
######################
# Loop to get consesus sequences
for R1_fastq, R2_fastq in zip(R1_fastqs, R2_fastqs):
        R1_barcode_list, R2_barcode_list = create_barcode_lists(R1_fastq, R2_fastq)

        # Make sequences into Seq objects
        R1_barcode_sequences = [Seq(barcode) for barcode in R1_barcode_list]
        R2_barcode_sequences = [Seq(barcode) for barcode in R2_barcode_list]

        # Make motif objects
        R1_motifs = motifs.create(R1_barcode_sequences)
        R2_motifs = motifs.create(R2_barcode_sequences)

        # Get consensus sequences
        R1_consensus = print(R1_motifs.consensus)
        R2_consensus = print(R2_motifs.consensus)

         # Get library IDs
        R1_id = R1_fastq.split('_')[0]
        R2_id = R2_fastq.split('_')[0]

        # Save consensus results
        results = pd.DataFrame([{'R1_library_ID':R1_id, 'R1_consensus':R1_consensus,
                                 'R2_library_ID':R2_id, 'R2_consensus': R2_consensus}])
        
        # Concat results to final dataframe
        consensus_df = pd.concat([consensus_df, results])

        # Make weblogo images
        R1_motifs.weblogo(f'{weblogo_path}/{R1_id}_R1_weblogo.png')
        R2_motifs.weblogo(f'{weblogo_path}/{R2_id}_R1_weblogo.png')
##################
# Save final results as csv
outfile_name = 'UBL_consensus_barcodes.csv'
outfile_save = input_dir+'/'+outfile_name
consensus_df.to_csv(outfile_save)

print(f'Consensus results are saved here:{outfile_save}')