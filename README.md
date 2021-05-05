# OV-Project-One
Contains the code and scripts used to process the data for my OVCAR3 analysis 
# CRISPRi Screen Analysis
## Pre-processing - Screen_Preprocessing.R
The script used to perform VST normalization and batch correction for the 96 well ovcar3 screen data is called Screen_Preprocessing.R. This file took the raw counts as input, identified two complete outlier samples (which were removed), then performed subsequent normalization and batch correction using surrogate variable analysis; the output was then used for rank based analysis. 
## Rank Based Analysis - OVCAR3_Screen_Rank_Based_Analysis.ipynb
The output from pre-processing was fed into the jupyter notebook OVCAR3_Screen_Rank_Based_Analysis.ipynb. As we only had single replicates this method, I developed a method loosly based off of the CMAP project (CITE) to determine potential regulated genes for my super-enhancers. In this script every gene was ranked based its mean expression across all samples and change in rank values were determine in a gene-sample specific manner. I then used an internal program, written to determine the eFDR of detected genes at a given change in rank threshold, to determine regulated genes at a given emprical FDR.  
