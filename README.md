# OV-Project-One
Contains the code and scripts used to process the data for my OVCAR3 analysis 
# CRISPRi Screen Analysis
## Pre-processing - Screen_Preprocessing.R
The script used to perform VST normalization and batch correction for the 96 well ovcar3 screen data is called Screen_Preprocessing.R. This file took the raw counts as input, identified two complete outlier samples (which were removed), then performed subsequent normalization and batch correction using surrogate variable analysis; the output was then used for rank based analysis. 
## Rank Based Analysis - OVCAR3_Screen_Rank_Based_Analysis.ipynb
The output from pre-processing was fed into the jupyter notebook OVCAR3_Screen_Rank_Based_Analysis.ipynb. As we only had single replicates this method, I developed a method loosly based off of the CMAP project (CITE) to determine potential regulated genes for my super-enhancers. In this script every gene was ranked based its mean expression across all samples and change in rank values were determine in a gene-sample specific manner. I then used an internal program, written to determine the eFDR of detected genes at a given change in rank threshold, to determine regulated genes at a given emprical FDR.  
# Cancer Copy Number Analysis 
## Assignment Genome Wide - OVLP_CNV_Whole_Genome.py
This program was written to assign copy number from TCGA data to uniform 15kb sliding windows across the genome. This is done one chromosome at a time due to the size and complexity of the data; the shell script launch_pybatch_genome.sh written for a CENTOS architexture allows for parallel analysis of all chromosomes at one time. Input files and output files need to be adjusted as needed. 
## Assignment to SE Regions - SEOVLP_CNV.py
This program was written to assign copy number from TCGA data to uniform 15kb sliding windows across the genome which overlap super-enhaners (as determined by bedtools intersect). This is done one chromosome at a time due to the size and complexity of the data; the shell script run_pybatch_SEOVLP.sh written for a CENTOS architexture allows for parallel analysis of all chromosomes at one time. Input files and output files need to be adjusted as needed. 
## Combination - Combine_CNV_Chr_Files.ipynb
There are many ways to combine the outputs from these chromosome assignments into a whole genome subset, this is just the approach I used. It is rough but gets the job done. 
## Comparisons
## Amplifcation of super-enhancers - OVCAR_CNV_Comparison_Final.R
This program takes the combined output of copy number across the genome from all patients as well as a list of super-enhancer overlapping regions and performs a comparison of the amplification of the super-enhancer windows vs the genome via random subsetting. This is done using lopping coupled with random window subsets and a one sided t.test comparing SE overlapping windows to the randomly drawn genomic windows.  
## Copy Number eQTL analysis - CNV_eQTL.R
Herein I make use of MatrixEQTL and permutation testing to determine the number of significant Copy Number Associated Expression Quantitative Trait Loci in my ovarian cancer super-enhancer data; using TCGA data. As MatrixEQTL was built to compare variants (qualitative) to gene expression, and I am comparing copy number (quantitative) to gene expression I used an empirical false discovery rate to determine "significance" rather than a bonferonni method. In this program first the true experimental data is analyzed for signficant findings; then the RNA data is permutated at the column level (columns are moved around) while the CNV data is untouched. MatrixEQTL is then run on the permutated "null" data n times to determine the eFDR.  
