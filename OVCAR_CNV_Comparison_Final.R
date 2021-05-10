library(Matrix)
library(ggplot2)
library(matrixStats)
library(data.table)


setwd('/Users/mkelly95/Documents/Scripts_for_Projects/CNV_Stuff/CNV_Comparison/')

#whole genome sliding windows
CNV_fullset <- fread('ovcar_regions_with_OV_values_chr1to22.tsv', sep = '\t', header = TRUE)
#regions list
regions<- CNV_fullset$chromosome_identifier
regions<-regions[1:192080]

#read in the SE regions  
CNVSE_OVLP <- fread('Full_Genome_CNV_15kb_Sliding_Windows_SEOVLP_Regions.txt', sep = '\t', header = FALSE)
#subset the regions in my se overlap
CNVOverlapset<-CNV_fullset[CNV_fullset$chromosome_identifier %in% CNVSE_OVLP$V4]

tCNV <- t(CNVOverlapset)
tCNVnolabs<- tCNV[6:572,]

subregsn<- sample(regions,336,replace = FALSE)
CNVOverlapset_randomN<-CNV_fullset[CNV_fullset$chromosome_identifier %in% subregsn]
tCNVRn<- t(CNVOverlapset_randomN)
tCNVRnnolabs<- tCNVRn[6:572,]
#unpaired one sided t.test looking for amplification of the SE regions vs randomly drawn set
ttestres<-t.test(as.numeric(tCNVnolabs),as.numeric(tCNVRnnolabs), paired = FALSE, alternative = 'greater')


CNV_regions<-CNVSE_OVLP$V4
split_sub_regs<- strsplit(CNV_regions,'_')
#initilize
start_stop<-as.numeric(unlist(split_sub_regs[1])[2:3])
length_of_regions1<-start_stop[2]-start_stop[1]
all_start_stop_lengths<- c(length_of_regions1)
for (i in c(2:336)){
  start_stop<-as.numeric(unlist(split_sub_regs[i])[2:3])
  length_of_regions<-start_stop[2]-start_stop[1]
  all_start_stop_lengths<- c(all_start_stop_lengths,length_of_regions)}
#median length
cnv_regions_size_median<-median(all_start_stop_lengths)





#initialize t.test matrices to fill 
ttestpvals<- c(ttestres$p.value)
ttestest<- c(ttestres$estimate)
xest<-ttestres$estimate[1]
yest<-ttestres$estimate[2]
#make a uniform probability matrix
probof1 <-1/192080
probarray<-rep(probof1,192080)
#iterate, set number to preferred runs 
for (i in c(1:100))
{
  #track
  print(i)
  subregsn<- sample(regions,336,replace = FALSE, prob = probarray)
  #median stuff
  split_sub_regs<- strsplit(subregsn,'_')
  #get lengths need 1st region to initialize matrix
  start_stop<-as.numeric(unlist(split_sub_regs[1])[2:3])
  length_of_regions1<-start_stop[2]-start_stop[1]
  all_start_stop_lengths<- c(length_of_regions1)
  #for the remaining spots
  for (j in c(2:336)){
    start_stop<-as.numeric(unlist(split_sub_regs[j])[2:3])
    length_of_regions<-start_stop[2]-start_stop[1]
    all_start_stop_lengths<- c(all_start_stop_lengths,length_of_regions)}
  #median length check if the conditions are not met skip that i and get a new one
  if (median(all_start_stop_lengths) < cnv_regions_size_median*.95 || median(all_start_stop_lengths) > cnv_regions_size_median*1.05){next}
  #if these regions pass the check get the full data and t test it
  CNVOverlapset_randomN<-CNV_fullset[CNV_fullset$chromosome_identifier %in% subregsn]
  tCNVRn<- t(CNVOverlapset_randomN)
  tCNVRnnolabs<- tCNVRn[6:572,]
  #unpaired one sided t.test looking for amplification of the SE regions vs randomly drawn set
  ttestres<-t.test(as.numeric(tCNVnolabs),as.numeric(tCNVRnnolabs), paired = FALSE, alternative = 'greater')
  
  
  #fill arrays 
  ttestpvals<- c(ttestpvals,ttestres$p.value)
  #unpaired estimates 
  xest<- c(xest,ttestres$estimate[1])
  yest<- c(yest,ttestres$estimate[2])
}

#xest is the SE region values, yest is the random sample values
#therefore x-y is the difference of SE vs random
ovcarest<-xest-yest

#interesting stats 
sum(ttestpvals)
length(ttestpvals)
#the amount of positive estimates is the number of cases where SE amplification > random subset
sum(ovcarest >0)
#how many p values are 0
sum(ttestpvals == 0)
#how many are significant 
sum(ttestpvals < 0.05)


#summary plots

plotstuffbox<-boxplot(as.numeric(tCNVnolabs),as.numeric(tCNVRnnolabs),names= c('SE CNVs','Non-SE CNVs Random Subset '),ylim=c(-1,1),
                      col=c('orchid','white'))

plotstuffbox<-boxplot(as.numeric(xest),as.numeric(yest),names= c('SE CNV Est. Means','Non-SE CNVs Est. Means '),ylim=c(-.1,.1),
                      col=c('orchid','white'))
