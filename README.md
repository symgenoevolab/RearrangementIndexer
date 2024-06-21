![RearrangementIndexer](https://github.com/symgenoevolab/RearrangementIndexer/assets/143068437/a22291b2-2c0c-41a3-bce7-52d0f5db9af1)

# RearrangementIndexer.py: Quantify ancestral linkage group rearrangement in chromosome-level genomes

This script takes as input a directory of 'coordinates files' containing information about genes' genomic location and ALG, and calculates a 'Rearrangement Index' for each file.

Associated manuscript: 
Lewin TD, Liao IJY, Luo YJ. Annelid comparative genomics and the evolution of massive lineage-specific genome rearrangement in bilaterians. BioRxiv (2024). (https://www.biorxiv.org/content/10.1101/2024.05.15.594353v1) Please cite this manuscript if you use this script. 

## Should I use this script?

If you have a set of bilaterian genomes and you want to know to what extent the genes within them have been scrambled by interchromosomal rearrangements, YES. For anything else, NO.

## Executing the script

```
python RearrangementIndexer.py input_directory
```

* This script was tested using python version 3.10.12. It requires the data manipulation tool pandas. 

## Detailed description

### Script rationale:
The ancestor of bilaterians had 24 chromosomes. Genes that were together on these chromosomes are known as bilaterian ‘ancestral linkage groups’ (ALGs): these were inherited by all bilaterians (Simakov et al. 2022). Genome rearrangements such as chromosome fusion and fission may combine or split ALGs. Despite this, species from phyla as diverse as chordates molluscs and annelids have maintained this conserved genome structure for over half a billion years since the ancestor of bilaterians (Simakov et al 2013). This suggests that it has significant importance to genome function.

We developed a Rearrangement Index (Ri) to quantify the extent of ALG reaarrangement in bilaterian genomes. It is an very simple calculation that provides a summary of (a) whether genes from each ALG are kept together or split apart (a proxy for fission) and (b) whether they are combined with genes from other ALGs or kept in isolation (a proxy for fusion).

### Calculation of index:
For each ALG, the rearrangement index is calculated as:

                 RALG = 1 - (SCHR x CCHR)    

where RALG denotes the rearrangement index for a given ALG; 
SCHR (ALG splitting parameter) represents the highest proportion of genes from this ALG on a single chromosome;
CCHR (ALG combining parameter) is the proportion of genes on that chromosome that belong to that particular ALG. 

By incorporating these parameters, the index accounts for both ALG splitting and ALG combining. 

Subsequently, the Ri for each genome is given by the equation:
                   
                Ri =  (SUM(RALG))/N

where Ri denotes the rearrangement index for the genome; 
RALG the is rearrangement index for each ALG; 
and N is the total number of ALGs. 

The higher the index, the higher the level of interchromosomal rearrangements.  


### Explanatory figure:
![Figure_S14_Supplementary_Ri_explanation_with_tables](https://github.com/symgenoevolab/RearrangementIndexer/assets/143068437/2fe49d8e-c372-41b2-abd2-bd295a33f739)
Explanatory figure | Explanation of the rearrangement index. Parts (A) to (F) each show a highly simplified example of how the rearrangement index is calculated under different amounts of ALG splitting and ALG combining. Grey horizontal lines represent chromosomes and colored vertical lines represent genes. In this scenario there are four ALGs (blue, green, yellow, red) and each has 5 genes. This does not vary in any of the examples.

A high rearrangement index indicates many interchromosomal rearrangements while a low index indicates few rearrangements. 

For instance, the RALG for example A (no rearrangements) blue ALG is calculated as follows: 5/5 blue genes are on one chromosome, so SCHR = 1; then, on this chromosome, 5/5 genes are blue genes so CCHR = 1. RALG = 1 – (1 × 1) = 0. In example A, this is true for all four ALGs. So Ri = (0 + 0 + 0 + 0) / 4 = 0. 

In example B (ALG splitting by chromosome fission), the chromosome containing the blue ALG has fissioned into two parts. In this case, 3/5 blue genes are on one chromosome, so SCHR = 0.6; then, on this chromosome, 3/3 genes are blue genes (because there has been no fusion) so CCHR = 1. RALG = 1 – (0.6 × 1) = 0.4. In example B, the other three ALGs have RALG = 0, so Ri = (0.4 + 0 + 0 + 0) / 4 = 0.1.  

In example D (ALG combining by chromosome fusion), the chromosomes containing the blue and green ALGs have fused. For all ALGs SCHR = 1 because all genes are on the same chromosome. But, CCHR = 0.5 for the green and blue ALGs because only 5/10 of the genes on this new chromosome are from each one. Therefore, for green and blue ALGs, RALG = 1 – (0.5 × 1) = 0.5. For red and yellow ALGs, RALG = 1 – (1 × 1) = 0. Then, Ri = (0.5 + 0.5 + 0 + 0) / 4 = 0.25.  

Example G is the most complicated scenario, including extensive ALG splitting and combining similar to that observed in clitellate annelids. We will take the ALGs one by one. For the blue ALG, 3/5 genes are on one chromosome, so SCHR = 0.6. On this chromosome, 3/4 genes are blue, so CCHR = 0.75. RALG = 1 – (0.6 × 0.75) = 0.55. For the green ALG, 2/5 genes are on one chromosome, so SCHR = 0.4. On this chromosome, 2/4 genes are green, so CCHR = 0.5. RALG = 1 – (0.4 × 0.5) = 0.8. For the yellow ALG, 2/5 genes are on one chromosome, so SCHR = 0.4. On this chromosome, 2/4 genes are yellow, so CCHR = 0.5. RALG = 1 – (0.4 × 0.5) = 0.8. Finally, for the red ALG, 2/5 genes are on one chromosome, so SCHR = 0.4. On this chromosome, 2/2 genes are red, so CCHR = 0.4. RALG = 1 – (0.4 × 1) = 0.6. Then, Ri = (0.55 + 0.8 + 0.8 + 0.6) / 4 = 0.69.  This represents the highest rearrangement index in this simple set of scenarios. 

### Input files:

The script expects a directory containing tsv files with the following format:

```
1       Complete        CAIIXF020000008.1       36937317        36938462        A1b
2       Complete        CAIIXF020000005.1       32633780        32633985        Eb
3       Complete        CAIIXF020000001.1       61230213        61230325        P
4       Complete        CAIIXF020000007.1       27374167        27374209        B3
5       Complete        CAIIXF020000009.1       30922059        30922235        H
6       Complete        CAIIXF020000012.1       11914594        11914622        G
```

where each row is one gene and:
column 1 = Gene ID 
column 2 = Status of gene (not important for rest of script)
column 3 = Chromosome ID
column 4 = Gene start
column 5 = Gene end
column 6 = ALG

Genes for each species should be placed in separate tsv files.

How can I produce these files?  The output of our macrosynteny pipeline SyntenyFinder (https://github.com/symgenoevolab/SyntenyFinder) produces files called Species_coordinates.tsv which fit this format. However, as long as the files include the appropriate information, the method of generation does not matter.

### Output files:

Rearrangement_index.tsv = Rearrangement index for each ALG in each species.

Splitting_parameter.tsv = Splitting parameter for each ALG in each species (Splitting index = 1 splitting parameter).

Combining_parameter.tsv = Combining parameter for each ALG in each species (Combining index = 1 combining parameter).

## Version History

* 0.2
    * Various bug fixes and optimizations
* 0.1
    * Initial Release
