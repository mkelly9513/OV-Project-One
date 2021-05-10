library(Matrix)
library(ggplot2)
library(matrixStats)
library(MatrixEQTL)
library(data.table)
library(permute)
library(picante)


setwd('/Users/mkelly95/Documents/Scripts_for_Projects/CNV_Stuff/eQTL')
#loading base data 
useModel = modelLINEAR;
#CNV values 
SNP_file_name = paste( "OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_with_rnaseq_cnvfloatvals.txt", sep="'\t");
#Cleaned gene expression, removed genes with 100 0's or more 
expression_file_name = paste( "OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_RNA_SEQ_Data_less100NAns.tsv", sep="'\t");

#covariates_file_name = paste( "/data/Covariates.txt", sep="");
output_file_name = tempfile();


#loose
pvOutputThreshold = 1e-3

errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

#RUN EQTL
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  #cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);


## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
True_eQTLs<-dim(me$all$eqtls)[1]

plot(hist(me$all$eqtls$pvalue,breaks = 10))





#permutation of NULL
#adjust the number in seq to change the number of loops 
totalQTLCount <- c()
for (i in seq(1,10)) {
  dttest <- fread('OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_RNA_SEQ_Data_less100NAns.tsv', sep = '\t')
  df4 <- dttest
  mixedset <- shuffle(ncol(df4)-1)
  mixed1 <- mixedset+1
  newindc <- c(1)
  newindc <- append(newindc,mixed1)
  df4<-df4[,..newindc]
  write.table(df4,file='OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_RNA_SEQ_col_shuffle_move.tsv',sep = '\t', row.names = FALSE)
  
  
  expression_file_name = paste("OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_RNA_SEQ_col_shuffle_move.tsv", sep = "'\t");
  SNP_file_name = paste( "OVCAR3_Super_Enhancer_Overlap_15kbwindows_patients_with_rnaseq_cnvfloatvals.txt", sep="'\t");
  
  #pvOutputThreshold = 1.594388e-11;
  
  
  
  
  #loose
  pvOutputThreshold = 1e-3
  
  errorCovariance = numeric();
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    #cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name);
  
  #cor(dttest[,2:292],)
  ## Results:
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected eQTLs:', '\n');
  
  numeQTLs<-dim(me$all$eqtls)[1]
  totalQTLCount <- append(totalQTLCount,numeQTLs)
  
}

plot(hist(me$all$eqtls$pvalue,breaks = 10))

#Significance
#we take the median null CNeQTL to serve
mednull<- median(totalQTLCount)
#Then divide this by the true number of CNeQTLs to determine the empirical FDR
FDRval<- mednull/True_eQTLs



