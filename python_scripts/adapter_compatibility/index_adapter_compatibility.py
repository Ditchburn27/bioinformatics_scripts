import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import textwrap

# Function to calculate the Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

def main():
    # Set up command-line argument parsing
    description = textwrap.dedent('''\
        Generate a Hamming distance heatmap from two CSV files,
        containing index names in 1st column & index sequences in 2nd column.

        Version: 1.0
        by Leighton Ditchburn
        ''')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f1', '--file1', type=str, required=True, help='Path to the first CSV file')
    parser.add_argument('-f2', '--file2', type=str, required=True, help='Path to the second CSV file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Name of the output file for the heatmap')
    parser.add_argument('-t', '--title', type=str, default='Hamming Distance Heatmap', help='Title of the heatmap plot')
    args = parser.parse_args()
    
    # Extract base names for labels
    base_name1 = os.path.basename(args.file1).split('.')[0]
    base_name2 = os.path.basename(args.file2).split('.')[0]

    # Read the two CSV files
    df1 = pd.read_csv(args.file1, header=0, dtype={0: str, 1: str})  # Read first row as header, ensure columns are strings
    df2 = pd.read_csv(args.file2, header=0, dtype={0: str, 1: str})  # Read first row as header, ensure columns are strings

    # Remove duplicate rows
    df1 = df1.drop_duplicates()
    df2 = df2.drop_duplicates()

    # Ensure the columns are named appropriately
    name_col_name = df1.columns[0]  # First column as name
    index_col_name = df1.columns[1]  # Second column as index

    names1 = df1[name_col_name].tolist()
    names2 = df2[name_col_name].tolist()
    indexes1 = df1[index_col_name].tolist()
    indexes2 = df2[index_col_name].tolist()

    # Clean up index sequences (remove leading/trailing whitespace)
    indexes1 = [idx.strip() if isinstance(idx, str) else '' for idx in indexes1]
    indexes2 = [idx.strip() if isinstance(idx, str) else '' for idx in indexes2]

    # Print out the sequences for debugging
    print("Indexes1:", indexes1)
    print("Indexes2:", indexes2)

    # Check for empty sequences and print them out
    for idx in indexes1:
        if not idx:
            print("Empty sequence found in file1:", idx)
    for idx in indexes2:
        if not idx:
            print("Empty sequence found in file2:", idx)

    # Calculate the Hamming distance matrix
    hamming_matrix = np.zeros((len(indexes1), len(indexes2)))

    for i, idx1 in enumerate(indexes1):
        for j, idx2 in enumerate(indexes2):
            try:
                hamming_matrix[i, j] = hamming_distance(idx1, idx2)
            except ValueError as e:
                print(f"Error calculating Hamming distance between {idx1} and {idx2}: {e}")
                return

    # Filter out rows and columns without any Hamming distances less than 5
    row_mask = np.any(hamming_matrix < 5, axis=1)
    col_mask = np.any(hamming_matrix < 5, axis=0)
    filtered_hamming_matrix = hamming_matrix[row_mask][:, col_mask]

    filtered_names1 = [f"{n}-{i}" for n, i, keep in zip(names1, indexes1, row_mask) if keep]
    filtered_names2 = [f"{n}-{i}" for n, i, keep in zip(names2, indexes2, col_mask) if keep]

    # Create a DataFrame for the heatmap
    hamming_df = pd.DataFrame(filtered_hamming_matrix, index=filtered_names1, columns=filtered_names2)

    # Plot the heatmap with red-green color scale
    plt.figure(figsize=(30, 24))
    sns.heatmap(hamming_df, annot=True, cmap='RdYlGn', cbar=True, fmt='.0f')
    plt.title(args.title)  # Use the provided title for the plot
    plt.xlabel(base_name2)  # Use base names for x-axis
    plt.ylabel(base_name1)  # Use base names for y-axis
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    
    # Adjust margins to make more space for the labels
    plt.subplots_adjust(bottom=0.3, left=0.3)

    plt.savefig(f'{args.output}.pdf')  # Save the heatmap to the specified output file
    plt.show()

if __name__ == '__main__':
    main()
