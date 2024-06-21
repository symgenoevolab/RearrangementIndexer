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

### Subsequently, the Ri for each genome is given by the equation:
                      
                Ri =  (SUM(RALG))/N

where Ri denotes the rearrangement index for the genome; 
RALG the is rearrangement index for each ALG; 
and N is the total number of ALGs. 

The higher the index, the higher the level of interchromosomal rearrangements.  

### Input files

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

### Output files 
Rearrangement_index.tsv = Rearrangement index for each ALG in each species.

Splitting_parameter.tsv = Splitting parameter for each ALG in each species (Splitting index = 1 splitting parameter).

Combining_parameter.tsv = Combining parameter for each ALG in each species (Combining index = 1 combining parameter).

## Version History

* 0.2
    * Various bug fixes and optimizations
* 0.1
    * Initial Release
