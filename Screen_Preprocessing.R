library(Matrix)
library(ggplot2)
library(matrixStats)
library(data.table)
library(sva)
library(DESeq2)
library(ggplot2)

####Basics
#write.table(OVCAR3lowzero.dt,file='OVCAR3_96_Well_Screen_LowExpressed_Genes_Removed.tsv',sep = '\t', row.names = FALSE)

#design matrix
Design.dt<- fread('DesignMatrix.txt', sep = '\t', header = TRUE)
Design.dt$Day <- as.factor(Design.dt$Day)
Design.dt$Condition <- as.factor(Design.dt$Condition)
Design.dt$Block <- as.factor(Design.dt$Block)
#DESEQ2

########

##############
# raw counts
DESOV<- fread('OVCAR3_96_Well_Screen_LowExpressed_Genes_Removed.tsv', sep = '\t', header = TRUE)
OVCAR3lowzero.dt <-fread('OVCAR3_96_Well_Screen_LowExpressed_Genes_Removed.tsv', sep = '\t', header = TRUE)
rownames(DESOV) <- c(DESOV$Gene)
DESOV$Gene <- NULL
dds <- DESeqDataSetFromMatrix(countData = DESOV,
                              colData = Design.dt,
                              design= ~ Condition + Day +Block )

dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
svd <- vst(dds, blind = FALSE,nsub = 4000)

pcaData <- plotPCA(svd, intgroup=c("Condition", "Block"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Block)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv))   n.sv <- num.sv(dat,mmi,method="leek")
  o <- svaseq(dat,mmi,mm0,n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)}


#setup need to adjust moba/b for each datap
datap <- assay(svd)
#datap<- assay(dds6164)

#look into non-vst corrected 
idx <- rowMeans(datap) > 1
dat <- datap[idx,]
moda <- model.matrix(~ Condition + Day, colData(svd))
modb <- model.matrix(~ 1, colData(svd))
correcteddata<- svaBatchCor(datap,moda,modb,3)
setest<- SummarizedExperiment(correcteddata$corrected, colData = Design.dt)
correctedse<- DESeqTransform(setest)
rownames(correctedse) <- OVCAR3lowzero.dt$Gene


pcaData <- plotPCA(correctedse,  intgroup=c("Condition", "Day"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#Final Approach
#####################################
#####################################
#outliers are 61 and 64 and 62
#removed below
Design6164.dt<- fread('DesignMatrix_no6164.txt', sep = '\t', header = TRUE)
Design6164.dt$Day <- as.factor(Design6164.dt$Day)
Design6164.dt$Condition <- as.factor(Design6164.dt$Condition)
Design6164.dt$Block <- as.factor(Design6164.dt$Block)
Design6164.dt$Plate <- as.factor(Design6164.dt$Plate)


DESOV61<- OVCAR3lowzero.dt
rownames(DESOV61) <- c(DESOV61$Gene)
DESOV61$Gene <- NULL
DESOV61$OVCAR3_H7_61 <- NULL
DESOV61$OVCAR3_F2_64 <- NULL
dds6164 <- DESeqDataSetFromMatrix(countData = DESOV61,
                              colData = Design6164.dt,
                              design= ~ Condition + Day +Block )

dds6164<-estimateSizeFactors(dds6164)
dds6164<-estimateDispersions(dds6164)
svd6164 <- vst(object = dds6164, blind = FALSE,nsub = 5000)
pcaData6164 <- plotPCA(svd6164, intgroup=c("Condition", "Day"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData6164, "percentVar"))
ggplot(pcaData6164, aes(PC1, PC2, color=Condition, shape=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
#################remove SVA stuff
svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv))   n.sv <- num.sv(dat,mmi,method="leek")
  o <- svaseq(dat,mmi,mm0,n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)}


#setup need to adjust moba/b for each datap
datap <- assay(svd6164)
#datap<- assay(dds6164)

#diffmodel
moda <- model.matrix(~ Block + Day , colData(svd6164))


modb <- model.matrix(~ 1, colData(svd6164))


correcteddata<- svaBatchCor(datap,moda,modb,2)


setest<- SummarizedExperiment(correcteddata$corrected, colData = Design6164.dt)
correctedse<- DESeqTransform(setest)
rownames(correctedse) <- OVCAR3lowzero.dt$Gene


pcaData <- plotPCA(correctedse,  intgroup=c("Condition", "Day"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition, shape=Day)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

#take corrected data 
testassay<- assay(correctedse)

testdf<- data.table(testassay)
testdf$Gene <- OVCAR3lowzero.dt$Gene
rownames(testdf) <-testdf$Gene

pythondt<- fread('OVCAR3_DT_With_Medians_Chromosomes_andstuff_for_loops.tsv', sep = '\t', header = TRUE)

DESEQ_Output_DT_noZero_count_Genes <- testdf[testdf$Gene %in% pythondt$Gene]

#write.table(DESEQ_Output_DT_noZero_count_Genes,file='DESEQ_Normalized_Output_dt_nonzero_genes.tsv',sep = '\t', row.names = FALSE)


