'''
This post processing script takes the output of the exon counting, merges the PCR batches and the random hexamer and shortdT reads from the same cells, and formats and saves the final output files.
'''
import sys
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.io as sio

## This function merges the different PCR batches, formats and saves the files
def merge_gene_count_files(input_folder, output_folder, sampleID, RT_matching_file, pipeline):
    
    ## Open sampleID file
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    ##Open gene annotation file
    gene_annotation_file = input_folder + "gene_name_annotate.csv"
    gene_annotation = pd.read_csv(gene_annotation_file, header=None)
    
    ## Open and merge the PCR batches
    for i, sample in enumerate(sample_list):
        cell_annotation_file = input_folder + sample + "_cell_annotate.csv"
        cell_annotation = pd.read_csv(cell_annotation_file, header=None)
        count_matrix_file = input_folder + sample + ".count"

        try:
            count_matrix = np.genfromtxt(count_matrix_file, delimiter=',')

            # Handle empty file - in case of 0 alignments
            if count_matrix.size == 0:
                print(f"Warning: {sample}.count is empty. Skipping.")
                continue

            # Reshape if only one line
            if count_matrix.ndim == 1:
                count_matrix = count_matrix.reshape(1, -1)

            sparse_matrix = sp.coo_matrix(
                (count_matrix[:, 2], ((count_matrix[:, 0] - 1), (count_matrix[:, 1] - 1))),
                shape=(gene_annotation.shape[0], cell_annotation.shape[0]),
                dtype=float
            )

            if i == 0 or 'final_count_matrix' not in locals():
                final_cell_annotation = cell_annotation
                final_count_matrix = sparse_matrix
            else:
                final_cell_annotation = pd.concat([final_cell_annotation, cell_annotation])
                final_count_matrix = sp.hstack([final_count_matrix, sparse_matrix])

        except Exception as e:
            print(f"Error processing {sample}: {e}. Skipping.")
            continue
    
    ## Divide the values by 2 if it is the paired-end pipeline, because this pipeline counts the 2 read pairs separately during exon counting
    if pipeline == "paired-end":
        final_count_matrix_corrected = final_count_matrix.tocsr() / 2
    else:
        final_count_matrix_corrected = final_count_matrix.tocsr()
    
    ## Add metadata to the cell annotation matrix and remove cells with 0 reads
    final_cell_annotation.columns = ['Cell_name', 'Unmatched_rate']
    final_cell_annotation['UMI_count'] = final_count_matrix_corrected.sum(axis=0).A[0]
    final_cell_annotation['Gene_count'] = (final_count_matrix_corrected > 0).sum(axis=0).A[0]
    final_cell_annotation[['PCR_batch', 'RT_Ligation_barcodes']] = final_cell_annotation['Cell_name'].str.split(pat = '.', n = 1, expand = True)
    final_cell_annotation['Ligation_barcodes'] = final_cell_annotation['RT_Ligation_barcodes'].str[0:10]
    final_cell_annotation['RT_barcodes'] = final_cell_annotation['RT_Ligation_barcodes'].str[10:20]
    final_count_matrix_corrected = final_count_matrix_corrected[:, np.where(final_cell_annotation['UMI_count'] > 0)[0]] 
    final_cell_annotation = final_cell_annotation[final_cell_annotation['UMI_count'] > 0].reset_index(drop=True)
    
    ## Format gene annotation dataframe
    final_gene_annotation = gene_annotation.iloc[:,[0,1,2,3]]
    final_gene_annotation.columns = ['exon_id', 'gene_type', 'gene_name', 'exon_number']
    
    ## Merge RT primers from one cell
    final_cell_annotation_merged, final_count_matrix_corrected_merged = merge_RT_primers(final_cell_annotation, final_count_matrix_corrected, RT_matching_file)
    
    ## Save files
    final_cell_annotation_merged.to_csv((output_folder + '/cell_annotation.csv'), index=False)
    final_gene_annotation.to_csv((output_folder + '/exon_annotation.csv'), index=False) 
    sio.mmwrite((output_folder + '/expression_matrix.mtx'), final_count_matrix_corrected_merged)

    
## This function merges the random hexamer and shortdT reads from the same cells.
def merge_RT_primers(cell_annotation, count_matrix, RT_matching_file):
    
    matching = pd.read_table(RT_matching_file, sep="\t", header=0)
    cell_annotation['RT_barcodes_new'] = cell_annotation['RT_barcodes']

    # Iterate over the rows of matching and replace values in cell_annotation['RT_barcodes_new'] to unify RT barcodes
    for i in range(len(matching)):
        index = cell_annotation['RT_barcodes_new'] == matching.iloc[i, 1]
        cell_annotation.loc[index, 'RT_barcodes_new'] = matching.iloc[i, 0]

    # Create 'New_cell_name' column with the unifies RT barcode across shortdT and randomN RT primers
    cell_annotation['New_cell_name'] = cell_annotation['PCR_batch'] + '.' + cell_annotation['Ligation_barcodes'] + cell_annotation['RT_barcodes_new']

    # Collapse same cells
    table_counts = cell_annotation['New_cell_name'].value_counts()
    singlet = table_counts[table_counts == 1].index
    doublet = table_counts[table_counts == 2].index

    # Separate singlet and doublet cells, then collapse doublet cells    
    cell_annotation_singlet = cell_annotation[cell_annotation['New_cell_name'].isin(singlet)].reset_index(drop=True)
    count_matrix_singlet = count_matrix[:, np.where(cell_annotation['New_cell_name'].isin(singlet))[0]]

    cell_annotation_singlet['ShortdT_UMI_count'] = 0
    cell_annotation_singlet['RandomN_UMI_count'] = 0
    cell_annotation_singlet['ShortdT_UMI_count'].loc[cell_annotation_singlet['RT_barcodes'].isin(matching.iloc[:, 0])] = cell_annotation_singlet['UMI_count'].loc[cell_annotation_singlet['RT_barcodes'].isin(matching.iloc[:, 0])]
    cell_annotation_singlet['RandomN_UMI_count'].loc[cell_annotation_singlet['RT_barcodes'].isin(matching.iloc[:, 1])] = cell_annotation_singlet['UMI_count'].loc[cell_annotation_singlet['RT_barcodes'].isin(matching.iloc[:, 1])]

    cell_annotation_doublet = cell_annotation[cell_annotation['New_cell_name'].isin(doublet)].reset_index(drop=True)
    cell_annotation_doublet_sorted = cell_annotation_doublet.sort_values('New_cell_name').reset_index(drop=True) 
    cell_annotation_doublet_first_cell = cell_annotation_doublet_sorted.iloc[::2].reset_index(drop=True)
    cell_annotation_doublet_second_cell = cell_annotation_doublet_sorted.iloc[1::2].reset_index(drop=True)

    cell_annotation_doublet_first_cell['ShortdT_UMI_count'] = 0
    cell_annotation_doublet_first_cell['RandomN_UMI_count'] = 0
    cell_annotation_doublet_first_cell['ShortdT_UMI_count'].loc[cell_annotation_doublet_first_cell['RT_barcodes'].isin(matching.iloc[:, 0])] = cell_annotation_doublet_first_cell['UMI_count'].loc[cell_annotation_doublet_first_cell['RT_barcodes'].isin(matching.iloc[:, 0])]
    cell_annotation_doublet_first_cell['RandomN_UMI_count'].loc[cell_annotation_doublet_first_cell['RT_barcodes'].isin(matching.iloc[:, 1])] = cell_annotation_doublet_first_cell['UMI_count'].loc[cell_annotation_doublet_first_cell['RT_barcodes'].isin(matching.iloc[:, 1])]
    cell_annotation_doublet_second_cell['ShortdT_UMI_count'] = 0
    cell_annotation_doublet_second_cell['RandomN_UMI_count'] = 0
    cell_annotation_doublet_second_cell['ShortdT_UMI_count'].loc[cell_annotation_doublet_second_cell['RT_barcodes'].isin(matching.iloc[:, 0])] = cell_annotation_doublet_second_cell['UMI_count'].loc[cell_annotation_doublet_second_cell['RT_barcodes'].isin(matching.iloc[:, 0])]
    cell_annotation_doublet_second_cell['RandomN_UMI_count'].loc[cell_annotation_doublet_second_cell['RT_barcodes'].isin(matching.iloc[:, 1])] = cell_annotation_doublet_second_cell['UMI_count'].loc[cell_annotation_doublet_second_cell['RT_barcodes'].isin(matching.iloc[:, 1])]
    cell_annotation_doublet_first_cell['ShortdT_UMI_count'] = cell_annotation_doublet_first_cell['ShortdT_UMI_count'] + cell_annotation_doublet_second_cell['ShortdT_UMI_count']
    cell_annotation_doublet_first_cell['RandomN_UMI_count'] = cell_annotation_doublet_first_cell['RandomN_UMI_count'] + cell_annotation_doublet_second_cell['RandomN_UMI_count']

    total_reads_first = 1 / (1 - cell_annotation_doublet_first_cell.iloc[:, 1]) * cell_annotation_doublet_first_cell.iloc[:, 2]
    total_reads_second  = 1 / (1 - cell_annotation_doublet_second_cell.iloc[:, 1]) * cell_annotation_doublet_second_cell.iloc[:, 2]
    cell_annotation_doublet_first_cell.iloc[:, 1] = (cell_annotation_doublet_first_cell.iloc[:, 1] * total_reads_first + cell_annotation_doublet_second_cell.iloc[:, 1] * total_reads_second) / (total_reads_first + total_reads_second)
    cell_annotation_doublet_first_cell.iloc[:, 2:4] = cell_annotation_doublet_first_cell.iloc[:, 2:4] + cell_annotation_doublet_second_cell.iloc[:, 2:4]
    count_matrix_doublet = count_matrix[:,np.where(cell_annotation['New_cell_name'].isin(doublet))[0]]   
    count_matrix_doublet_sorted = count_matrix_doublet[:, cell_annotation_doublet['New_cell_name'].sort_values().index]       
    count_matrix_doublet_first_cell = count_matrix_doublet_sorted[:, ::2]
    count_matrix_doublet_second_cell = count_matrix_doublet_sorted[:, 1::2]  
    count_matrix_doublet_merged = count_matrix_doublet_first_cell + count_matrix_doublet_second_cell
    cell_annotation_merged = pd.concat([cell_annotation_singlet, cell_annotation_doublet_first_cell])
    count_matrix_merged = sp.hstack([count_matrix_singlet, count_matrix_doublet_merged])
    cell_annotation_merged['Gene_count'] = (count_matrix_merged > 0).sum(axis=0).A[0]
    cell_annotation_merged['RT_barcodes_shortdT'] = cell_annotation_merged['RT_barcodes_new'] # RT barcode is consistently the shortdT barcode
    cell_annotation_merged['Cell_name'] = cell_annotation_merged['New_cell_name'] # Rename the cell to include consistently the shortdT barcode
    cell_annotation_merged = cell_annotation_merged.drop(columns=['RT_Ligation_barcodes', 'RT_barcodes_new', 'New_cell_name', 'RT_barcodes'])

    return cell_annotation_merged, count_matrix_merged
    


if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    sampleID = sys.argv[3]
    RT_matching_file = sys.argv[4]
    pipeline = sys.argv[5]
    merge_gene_count_files(input_folder, output_folder, sampleID, RT_matching_file, pipeline)


