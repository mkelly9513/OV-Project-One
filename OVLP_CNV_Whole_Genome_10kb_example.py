
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from argparse import ArgumentParser
#FYI, this program was written to analyze ovcar and BRCA data and the variables are named as such
#Ignore the naming of the variables and focus on what they do, they will be explained as they are assigned

# In[11]:


Parser = ArgumentParser(description='This will allow one to subset the input files for specific chromosome comparisons')
Parser.add_argument('--chrom', help='input "chr#" to subset regions file ',
                    type=str, required=True)
Parser.add_argument('--cnumber', help='input int number to get BRCA regions for that chromosome',
                    type=int, required=True)


# In[13]:


args = Parser.parse_args()


# In[12]:


ChromosomeForRegions= args.chrom
ChromNumberforBRCA= args.cnumber


# In[2]:

#testing
#ChromosomeForRegions= 'chr1'
#ChromNumberforBRCA= 1


# In[3]:


#insert regions of interest, ignore variable names, they apply to the original running of the data 
#regions file for whole genome
ovcar_regionsdt=pd.read_csv('Full_Genome_10kb_Sliding_Windows_hg19.txt', sep='\t', header=0)
#This is your data matrix of cnv regions, from TCGA
BRCA_data= pd.read_csv('OV_CNV_Data_Formatted.txt', header=0, sep='\t')



# In[9]:


#BRCA1['Sample'].iloc[15]


# In[4]:

#gather patient IDs
unique_ids= list(np.unique(BRCA_data['Sample']))


# In[5]:

#make your starting matrix to be replaced by true values, 0 aka diploid is assumed as the baseline
subzero = np.zeros(len(ovcar_regionsdt['start']))
for v in unique_ids:
    ovcar_regionsdt[v] = subzero


# In[6]:


ovcardtcolnames=list(ovcar_regionsdt.columns.values)


# In[7]:

#determine which samples are from tumors based on TCGA naming convention
tumorsamples = ['chromosome','start','stop']
acceptedset = ['01A','02A','03A','04A','05A','06A','07A','08A','09A']
subsettotest=ovcardtcolnames[4:]
for v in subsettotest:
    splitv = str.split(v,"-")
    if splitv[3] in acceptedset:
        tumorsamples.append(v)
ovcar_regionsdt_tumors = ovcar_regionsdt[tumorsamples]
ovcar_regionsdt_tumors=ovcar_regionsdt_tumors.reset_index(drop=True)

# In[15]:

#cut the data for a specific chromosome
BRCA1 = BRCA_data[BRCA_data['Chromosome'] == ChromNumberforBRCA]
BRCA1=BRCA1.reset_index(drop=True)
ovcardt1 = ovcar_regionsdt_tumors[ovcar_regionsdt_tumors['chromosome']== ChromosomeForRegions]
ovcardt1=ovcardt1.reset_index(drop=True) 
ovcardt1colnames=list(ovcardt1.columns.values)


# In[18]:


#ovcardt1



# In[12]:

#real loop that assigns signal to each region
#keep commented out lines commented out, they do not apply if using argument parser




for i in range(len(ovcardt1colnames)):
    if i > 2:
        colname= str(ovcardt1colnames[i])
        for j in range(len(BRCA1['Sample'])):
            #name = BRCA1['Sample'].iloc[j] 
            name = BRCA1.loc[j,'Sample'] 
            if name == colname:
                #copynumval = BRCA1['Segment_Mean'].iloc[j]
                copynumval = BRCA1.loc[j,'Segment_Mean']
                #copystart= BRCA1['Start'].iloc[j]
                copystart= BRCA1.loc[j,'Start']
                #copyend= BRCA1['End'].iloc[j]
                copyend= BRCA1.loc[j,'End']
                #chrom = BRCA_data['Chromosome'][j]
                for k in range(len(ovcardt1['start'])):
                    #boundstart= ovcardt1['start'].iloc[k]
                    boundstart= ovcardt1.loc[k,'start']
                    #boundend= ovcardt1['stop'].iloc[k]
                    boundend= ovcardt1.loc[k,'stop']
                    #chrom2 = int(str.split(ovcar_regionsdt['chromosome'][k],sep = 'chr')[1])
                    if boundstart >= copystart and boundend <= copyend:
                            #ovcardt1[colname].loc[k] = copynumval
                            ovcardt1.loc[k,colname] = copynumval




# In[ ]:

#output file
#alter the name to be what you want   
ovcardt1.to_csv(ChromosomeForRegions+'_OV_CNV_Tumor_Values_Over_Whole_Genome_Regions_10kb.csv', index=False)

