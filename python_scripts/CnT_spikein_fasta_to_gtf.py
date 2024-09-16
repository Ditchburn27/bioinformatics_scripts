# Script for converting SNAP-CUTANA spike-in FASTA file to GTF file
#######################################################
# import libraries
from Bio import SeqIO
import argparse

######################################################
def main():
    # Handle CLI arguments
    description = ('''\
            Convert FASTA reference file to GTF file with unique virtual chromosomes.

            Version: 1.1
            by Leighton Ditchburn
            ''')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output GTF file')
    args = parser.parse_args()

    fasta_file = args.input
    gtf_file = args.output

    #####################################################
    # Open the GTF file for writing
    with open(gtf_file, "w") as gtf:
        # Parse the FASTA file
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            # Create a unique virtual chromosome name based on the gene_id
            virtual_chromosome = f"chr{seq_record.id}"  # Assign virtual chromosome from gene_id
            
            source = "converted"  # The source of the annotation
            feature = "exon"  # Assume exon feature for simplicity
            start = 1  # Arbitrary start
            end = len(seq_record.seq)  # End based on the sequence length
            score = "."  # Placeholder for score
            strand = "+"  # Assume positive strand, modify as needed
            frame = "."  # Placeholder for frame
            
            # Extract identifiers for the attributes column from the FASTA header
            gene_id = seq_record.id  # Use FASTA header as the gene_id
            transcript_id = f"{seq_record.id}_t1"  # Create a transcript ID

            # Prepare the attributes field
            attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'

            # Write the GTF record with unique virtual chromosome
            gtf.write(f"{virtual_chromosome}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n")

    print(f"GTF file '{gtf_file}' created successfully with virtual chromosomes.")

#######################################################
if __name__ == '__main__':
    main()
