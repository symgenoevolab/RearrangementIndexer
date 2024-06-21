########################## RearrangementIndexer.py ###########################

### Author: Thomas D. Lewin
### Date created: 01/05/2024
### Last modified 20/06/2024
### Version: v0.2

############################ Brief Description ###############################

### Usage: python RearrangementIndexer.py input_directory

### Associated manuscript: 
### Lewin TD, Liao IJY, Luo YJ. Annelid comparative genomics and the evolution 
### of massive lineage-specific genome rearrangement in bilaterians. BioRxiv (2024)

### Summary: 
### This script takes as input a directory of 'coordinates files' containing 
### information about genes' genomic location and ALG, and calculates a 
### 'Rearrangement Index' for each file.

######################## Should I use this script? ##########################

### If you have a set of bilaterian genomes and you want to know to what 
### extent the genes within them have been scrambled by interchromosomal
### rearrangements, YES. For anything else, NO.

########################## Detailed Description #############################

### Script rationale:
### The ancestor of bilaterians had 24 chromosomes. Genes that were together 
### on these chromosomes are known as bilaterian ‘ancestral linkage groups’ 
### (ALGs): these were inherited by all bilaterians (Simakov et al. 2022).
### Genome rearrangements such as chromosome fusion and fission may combine 
### or split ALGs. Despite this, species from phyla as diverse as chordates,
### molluscs and annelids have maintained this conserved genome structure for
### over half a billion years since the ancestor of bilaterians (Simakov et
### al 2013). This suggests that it has significant importance to genome 
### function.
### We developed a Rearrangement Index (Ri) to quantify the extent of ALG
### reaarrangement in bilaterian genomes. It is an very simple calculation
### that provides a summary of (a) whether genes from each ALG are kept
### together or split apart and (b) whether they are combined with genes
### from other ALGs or kept in isolation.


### Calculation of index:
### For each ALG, the rearrangement index is calculated as:

                 # RALG = 1 - (SCHR x CCHR)    

### where RALG denotes the rearrangement index for a given ALG; SCHR 
### (ALG splitting parameter) represents the highest proportion of genes from
### this ALG on a single chromosome and CCHR (ALG combining parameter) is the
### proportion of genes on that chromosome that belong to that particular ALG. 
### By incorporating these parameters, the index accounts for both ALG 
### splitting and ALG combining. 

### Subsequently, the Ri for each genome is given by the equation:
                      
                # Ri =  (SUM(RALG))/N

### where Ri denotes the rearrangement index for the genome; 
### RALG the is rearrangement index for each ALG; 
### and N is the total number of ALGs. 
### The higher the index, the higher the level of interchromosomal
### rearrangements.  

############################## Input files #################################

### The script expects a directory containing tsv files with the following
### format:

### 1       Complete        CAIIXF020000008.1       36937317        36938462        A1b
### 2       Complete        CAIIXF020000005.1       32633780        32633985        Eb
### 3       Complete        CAIIXF020000001.1       61230213        61230325        P
### 4       Complete        CAIIXF020000007.1       27374167        27374209        B3
### 5       Complete        CAIIXF020000009.1       30922059        30922235        H
### 6       Complete        CAIIXF020000012.1       11914594        11914622        G

### where each row is one gene and:
### column 1 = Gene ID 
### column 2 = Status of gene (not important for rest of script)
### column 3 = Chromosome ID
### column 4 = Gene start
### column 5 = Gene end
### column 6 = ALG

### How can I produce these files? The output of our macrosynteny pipeline
### SyntenyFinder (https://github.com/symgenoevolab/SyntenyFinder) produces
### files called Species_coordinates.tsv which fit this format.

############################## Output files ################################

### Rearrangement_index.tsv = Rearrangement index for each ALG in each species.
### Splitting_parameter.tsv = Splitting parameter for each ALG in each species.
### Combining_parameter.tsv = Combining parameter for each ALG in each species.

############################# Start of script ##############################

# Import required packages
import pandas as pd
import os
import sys

# print citation message
print("If you use this script in your work, please cite: Lewin TD, Liao IJY, Luo YJ. Annelid comparative genomics and the evolution of massive lineage-specific genome rearrangement in bilaterians. BioRxiv (2024)")

# Define functions
def process_tsv_file(file_path):
    print(f"Processing file: {file_path}")
    
    # Load the TSV file into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=None)
    
    # Rename columns for clarity
    df.columns = ["Index", "Status", "Chromosome", "Start", "End", "ALG"]
    
    # Combine ALG sub-parts into single ALGs
    df['ALG'] = df['ALG'].replace({'A1a': 'A1', 'A1b': 'A1', 'Ea': "E", 'Eb': "E", 'Qa': 'Q', 'Qb': 'Q', 'Qc': 'Q', 'Qd': 'Q'})
    
    # Calculate proportion of genes belonging to each ALG for each chromosome
    chrom_alg_counts = df.groupby(["Chromosome", "ALG"]).size().unstack(fill_value=0)
    
    # Calculate the total count of genes for each chromosome
    chrom_total_counts = chrom_alg_counts.sum(axis=1)
    
    # Calculate the proportion of genes belonging to each ALG for each chromosome
    chrom_alg_proportions = chrom_alg_counts.div(chrom_total_counts, axis=0)
    
    # Determine all unique ALGs
    all_algs = chrom_alg_counts.columns
    
    # Determine chromosome with the most genes from each ALG
    chrom_most_genes = chrom_alg_counts.idxmax(axis=0)
    
    # Calculate the rearrangement index for each ALG
    rearrangement_indices = {}
    proportion_on_chrom_values = {}
    total_proportion_on_chrom_values = {}
    
    for alg in all_algs:
        most_genes_chrom = chrom_most_genes[alg]
        print(f"Calculating metrics for ALG: {alg}")
        
        # Find the proportion of genes from the ALG on the chromosome with the most genes from that ALG
        if alg in chrom_alg_counts.columns:
            proportion_on_chrom = chrom_alg_counts.loc[most_genes_chrom, alg] / chrom_alg_counts[alg].sum()
        else:
            proportion_on_chrom = pd.NA
        proportion_on_chrom_values[alg] = proportion_on_chrom
        
        # Find the proportion of total genes on that chromosome that are from that ALG
        if alg in chrom_alg_counts.columns:
            total_proportion_on_chrom = chrom_alg_counts.loc[most_genes_chrom, alg] / chrom_total_counts[most_genes_chrom]
        else:
            total_proportion_on_chrom = pd.NA
        total_proportion_on_chrom_values[alg] = total_proportion_on_chrom
        
        # Calculate the rearrangement index for this ALG
        if alg in chrom_alg_counts.columns:
            rearrangement_index = 1 - (proportion_on_chrom * total_proportion_on_chrom)
        else:
            rearrangement_index = pd.NA
        rearrangement_indices[alg] = rearrangement_index
    
    # Convert dictionaries to pandas Series
    rearrangement_series = pd.Series(rearrangement_indices, name=os.path.basename(file_path))
    proportion_on_chrom_series = pd.Series(proportion_on_chrom_values, name=os.path.basename(file_path))
    total_proportion_on_chrom_series = pd.Series(total_proportion_on_chrom_values, name=os.path.basename(file_path))
    
    return rearrangement_series, proportion_on_chrom_series, total_proportion_on_chrom_series

def main(input_dir):
    # Get list of TSV files in the input directory
    tsv_files = [file for file in os.listdir(input_dir) if file.endswith(".tsv")]
    
    # Initialize DataFrames to store the outputs
    rearrangement_df = pd.DataFrame()
    proportion_on_chrom_df = pd.DataFrame()
    total_proportion_on_chrom_df = pd.DataFrame()
    
    # Process each TSV file
    for tsv_file in tsv_files:
        tsv_path = os.path.join(input_dir, tsv_file)
        rearrangement_series, proportion_series, total_proportion_series = process_tsv_file(tsv_path)
        
        # Append to respective DataFrames
        rearrangement_df = pd.concat([rearrangement_df, rearrangement_series], axis=1)
        proportion_on_chrom_df = pd.concat([proportion_on_chrom_df, proportion_series], axis=1)
        total_proportion_on_chrom_df = pd.concat([total_proportion_on_chrom_df, total_proportion_series], axis=1)
    
    # Save the rearrangement index to a table format
    rearrangement_df.to_csv("Rearrangement_index.tsv", sep="\t")
    print("Rearrangement indices saved to Rearrangement_index.tsv")
    
    # Save the proportion_on_chrom values to a table format
    proportion_on_chrom_df.to_csv("Splitting_parameter.tsv", sep="\t")
    print("Fission parameters saved to Splitting_parameter.tsv")
    
    # Save the total_proportion_on_chrom values to a table format
    total_proportion_on_chrom_df.to_csv("Combining_parameter.tsv", sep="\t")
    print("Fusion parameters values saved to Combining_parameter.tsv")

if __name__ == "__main__":
    # Check if input directory is provided as command-line argument
    if len(sys.argv) != 2:
        print("Usage: python RearrangementIndexer.py input_directory")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    main(input_dir)

# print citation message
print("If you use this script in your work, please cite: Lewin TD, Liao IJY, Luo YJ. Annelid comparative genomics and the evolution of massive lineage-specific genome rearrangement in bilaterians. BioRxiv (2024)")

    ############################### End of script ##############################
