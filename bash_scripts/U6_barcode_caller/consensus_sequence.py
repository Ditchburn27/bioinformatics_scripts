import sys

print("Number of arguments:", len(sys.argv))
print("Argument List:", str(sys.argv))

def read_sequences(filename):
    with open(filename, 'r') as file:
        sequences = [line.strip() for line in file if line.strip()]
    return sequences

def consensus_sequence(sequences):
    # Assuming all sequences are of the same length
    length = len(sequences[0])
    consensus = []

    # Transpose the sequence list to process columns
    transposed = list(zip(*sequences))

    # For each position, find the most common base
    for column in transposed:
        consensus.append(max(set(column), key=column.count))

    return ''.join(consensus)

def update_reference(reference_file, consensus_sequence, barcode_type, library_id):
    with open(reference_file, 'r') as file:
        reference = file.read()

    # Determine the length of the random bases based on the barcode type
    if barcode_type == "borg" or barcode_type == "serloin":
        random_bases_length = 20
    elif barcode_type == "random":
        random_bases_length = 50
    else:
        print("Invalid barcode type")
        sys.exit(1)

    # Generate the string of Ns based on the length of random bases
    N_string = 'N' * random_bases_length

    # Replace the NNN bases with the consensus sequence
    updated_reference = reference.replace(N_string, consensus_sequence)

    # Construct the output filename with the library ID
    output_filename = f"{library_id}_{barcode_type}_reference.fa"

    with open(output_filename, 'w') as outfile:
        outfile.write(updated_reference)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_filename> <reference_filename> <barcode_type> <library_id>")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    reference_filename = sys.argv[2]
    barcode_type = sys.argv[3]
    library_id = sys.argv[4]

    sequences = read_sequences(input_filename)
    consensus = consensus_sequence(sequences)

    update_reference(reference_filename, consensus, barcode_type, library_id)