# OV-Project-One
Contains the code and scripts used to process and analyze the data for my analysis of super-enhancer function in ovarian cancer cells (OVCAR3) and in patient RNA-seq and Copy Number data from The Cancer Genome Atlas (TCGA). The OVCAR3 data for the CRISPR(i/KO) experiments can be downloaded from the GEO accession assocaited with our paper "A multi-omic dissection of super-enhancer driven oncogenic gene expression programs in ovarian cancer" ,(currently unpublished). Additionally there is OVCAR3 H3K9me3 and Hi-C Seq data. 
Single Cell Analysis was performed by Matt J Regner and can be found here:
Hi-C Seq Analysis was performed by Eric S Davis and can be found here:
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174259
# CRISPRi Screen Analysis
## Pre-processing - Screen_Preprocessing.R
The script used to perform VST normalization and batch correction for the 96 well ovcar3 screen data is called Screen_Preprocessing.R. This file took the raw counts as input, identified two complete outlier samples (which were removed), then performed subsequent normalization and batch correction using surrogate variable analysis; the output was then used for rank based analysis. 
## Rank Based Analysis - OVCAR3_Screen_Analysis_with_Plotting_LFC_Comparison.ipynb
### html: OVCAR3_Screen_Analysis_with_Plotting_LFC_Comparison_HTML.html
The output from pre-processing was fed into the jupyter notebook OVCAR3_Screen_Rank_Based_Analysis.ipynb. As we only had single replicates this method, I developed a method loosely based off of the CMAP project to determine potential regulated genes for my super-enhancers. In this script every gene was ranked based its mean expression across all samples and change in rank values were determine in a gene-sample specific manner. I then used an internal program, written to determine the eFDR of detected genes at a given change in rank threshold, to determine regulated genes at a given emprical FDR. I iterate through a number of changes in rank at a given FDR threshold to maximize detected DEGS. This code also includes various plotting functions and a log2 fold change comparitive analysis. There is an assocaited .hmtl readout (for those who just want to see the code and results without running it) called OVCAR3_Screen_Analysis_with_Plotting_LFC_Comparison_HTML.html . 
# CRISPR Knock Out Analysis 
## Pre-processing and DEA with DESEQ2 - DESEQ2_2021Reps_RNA_SVA_Plotting_V2.Rmd
### html: DESEQ2_2021Reps_RNA_SVA_Plotting_V2.html
In this R markdown file I walk through how I took the raw counts data (using a smaller subset toy example), built my input data matrix (full data at this point), cleaned the data, performed exploratory analysis, checked for batch correction, and then performed DEA. I used DESEQ2 to determine the genes that were differentially expressed (in this example just those that are downregulated) between wild type controls and SE14 knock out samples. I did this with a model accounting for batch internally in DESEQ2, with a base model (caring only about the KO vs WT), and a model accoutnign for one SV included. Any analysis performed here can be repeated for SE60 using the LFC and significance thresholds laid out in the corresponding paper. The raw counts data used as input for this analysis and the paper can be found at GEO. 

# Cancer Copy Number Analysis 
## Creating HG19 15kb windows - Split_Genome_into_windows.ipynb & Split_Genome_into_windows_10kb_example.ipynb
This program was used to take the hg19 reference genome chromosome lengths and create 15kb sliding windows for chromosomes 1-22 (autosomes). There is now an additional version called Split_Genome_into_windows_10kb_example.ipynb where I do this with 10kb windows. 
## Assignment Genome Wide - OVLP_CNV_Whole_Genome.py & OVLP_CNV_Whole_Genome_10kb_example.py
This program was written to assign copy number from TCGA data to uniform 15kb sliding windows across the genome. This is done one chromosome at a time due to the size and complexity of the data; the shell script launch_pybatch_genome.sh written for a CENTOS architexture allows for parallel analysis of all chromosomes at one time. Input files and output files need to be adjusted as needed. There is a new verison using the 10kb sliding windows as an example, depending on your version of python the syntax for assignment changes and the 10kb method has an updated assignment algorithm. 
## Assignment to SE Regions - SEOVLP_CNV.py
This program was written to assign copy number from TCGA data to uniform 15kb sliding windows across the genome which overlap super-enhaners (as determined by bedtools intersect). This is done one chromosome at a time due to the size and complexity of the data; the shell script run_pybatch_SEOVLP.sh written for a CENTOS architexture allows for parallel analysis of all chromosomes at one time. Input files and output files need to be adjusted as needed. 
## Combination - Combine_CNV_Chr_Files.ipynb & Combine_CNV_Chr_Files_10kb_example.ipynb
There are many ways to combine the outputs from these chromosome assignments into a whole genome subset, this is just the approach I used. It is rough but gets the job done. I added a new 10kb sliding window version to match the other file additions. 
## Comparisons
## Amplifcation of super-enhancers - OVCAR_CNV_Comparison_Final.R
This program takes the combined output of copy number across the genome from all patients as well as a list of super-enhancer overlapping regions and performs a comparison of the amplification of the super-enhancer windows vs the genome via random subsetting. This is done using lopping coupled with random window subsets and a one sided t.test comparing SE overlapping windows to the randomly drawn genomic windows.  
## Copy Number eQTL analysis - CNV_eQTL.R
Herein I make use of MatrixEQTL and permutation testing to determine the number of significant Copy Number Associated Expression Quantitative Trait Loci in my ovarian cancer super-enhancer data; using TCGA data. As MatrixEQTL was built to compare variants (qualitative) to gene expression, and I am comparing copy number (quantitative) to gene expression I used an empirical false discovery rate to determine "significance" rather than a bonferonni method. In this program first the true experimental data is analyzed for signficant findings; then the RNA data is permutated at the column level (columns are moved around) while the CNV data is untouched. MatrixEQTL is then run on the permutated "null" data n times to determine the eFDR. 
## Copy Number Survival Analysis - CNV_KM_Plots.R 
This program, written with the assistance of Dr. Joel Parker and Dr. Cheng "Chris" Fan, who created the internal function, takes survival data and copy number data in order to create KM plots and determine the hazard of amplification/deletion events across the genome. This was done with the curated TCGA survival data for all patients that had CNV data. 
