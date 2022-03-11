BiocManager::install("maftools")
library(maftools)
BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query)
sum_exp <- GDCprepare(query)
counts = assays(sum_exp)$"HTSeq - Counts"
counts
assays(sum_exp)$"HTSeq - Counts" [1:5,1:5]
str(sum_exp)
rowData(sum_exp)
head(rowData(sum_exp))
colData(sum_exp)
colData(sum_exp)[1:5, 25:29] #the rows are individual patients, while the columns are additional information about the patients
metadata(sum_exp)

#2.3
 dim(colData(sum_exp))
 dim(rowData(sum_exp))
 dim(assays(sum_exp)$"HTSeq - Counts")
 
#2.4
 str(colData(sum_exp))
 head(colData(sum_exp))
 #rows are the patient IDs and the columns are clinical information
 
#2.5
 colnames(colData(sum_exp))
 #age_at_diagnosis
 #age_at_index

 #2.6
 colData(sum_exp)$age_at_diagnosis[1:10]

 #2.7
 colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis / 365

 #2.8 
 colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis < 50, "Young", "Old")

 #2.9
 head(rowData(sum_exp))
 dim(rowData(sum_exp))
 #each column in rowData() represents alternate gene names. each row represents Ensembl Gene ID's 

 #2.10
 "MSH2" %in% rowData(sum_exp)$external_gene_name #results in TRUE
 "MLH1" %in% rowData(sum_exp)$external_gene_name #results in TRUE

 #2.11
 assays(sum_exp)$"HTSeq - Counts"[20:25,20:25]
 assays(sum_exp)$"HTSeq - Counts"[30:35,30:35]

 #2.12
 geneA_id_mask = (rowData(sum_exp)$external_gene_name == "MSH2")
 sum(geneA_id_mask)
 ensembl1_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask] 
 ensembl1_geneA 
 
 geneB_id_mask = (rowData(sum_exp)$external_gene_name == "MLH1")
 sum(geneB_id_mask)
 ensembl1_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask] 
 ensembl1_geneB 
 #both gene names match

 #2.13
 #the Ensembl gene ID is a row in rowData(sum_exp), not in assays(sum_exp)$"HTSeq - Counts"

 #2.14
 min(assays(sum_exp)$"HTSeq - Counts"[ensembl1_geneA,]) #the gene name goes on the left of the comma because rows are genes
 max(assays(sum_exp)$"HTSeq - Counts"[ensembl1_geneA,])
 summary(assays(sum_exp)$"HTSeq - Counts"[ensembl1_geneB,]) 
#2.15
 data_for_geneA = assays(sum_exp)$"HTSeq - Counts"["ENSG00000095002",]
 data_for_geneB = assays(sum_exp)$"HTSeq - Counts"["ENSG00000076242",]
 plot(data_for_geneA,
      data_for_geneB,
      xlab = "Gene A", 
      ylab = "Gene B"
  ) #there appears to be a positive correlation between the two genes
#2.16
 sum(is.na(colData(sum_exp)$age_category)) #if the NAs were to be removed, there owuld be less elements that indexed

 #2.17
 bool_age_na = is.na(colData(sum_exp)$age_category)
 num_na = sum(bool_age_na)
 num_na
 age_cat_no_NAs =  colData(sum_exp)$age_category[!bool_age_na]
 #2.18
 length(age_cat_no_NAs)
 dim(colData(sum_exp))[1] - num_na == length(age_cat_no_NAs)
 
 #2.19
 dim(assays(sum_exp)$"HTSeq - Counts")[2] #521 patients
 #this is more than the number of patients in age_cat_no_NAs because the NAs are included
 
 #2.20
 identical( rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")  ) #TRUE
 gene_counts = assays(sum_exp)$"HTSeq - Counts"["ENSG00000095002",!bool_age_na]
 #2.21
 length(gene_counts) == length(age_cat_no_NAs) #TRUE  
 #gene_counts was created by using a boolean mask to isolate the columns without NA values, that corresponded with the row of the gene name A
 #2.22
 boxplot(gene_counts ~ age_cat_no_NAs, 
         xlab = "Age Category", 
         ylab = "Gene Counts") #the older patients have many more gene counts in the upper quartile than the younger patients.
 
 #3.1
  #1) This owuld be accessed by assays(sum_exp)$"HTSeq - Counts". The rows represent genes and the columns represent Patient IDs.
  #2) Gene data can be accessed in the row data frame of the sum_exp data. This dataframe consists of the rows of the counts data. This row dataframe contains more info about the gene names.
  #3) Patient data can be accessed in the column data frame of the sum_exp data. This dataframe contains more informatio about the patients. This column dataframe consists of the rows of the counts data.