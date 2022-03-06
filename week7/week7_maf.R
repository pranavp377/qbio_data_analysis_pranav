#1.2
BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)
library(ggplot2)

#1.3
clinic <- data.table::fread("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data/coad_clinical_data.csv",
                            data.table = F)

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
colnames(clinic)
length(colnames(clinic)) #78
length(colnames(clinic) == "bcr_patient_barcode") #78, data type boolean
#ther are no TRUEs because the column name got changed

#1.4
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#1.5
getwd()
setwd("~/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data")
list.files()
setwd("~/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data/GDCdata")
list.files()
maf_dataframe <- data.table::fread("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                  data.table = F)
clinic <- data.table::fread("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data/coad_clinical_data.csv",
                            data.table = F)
#2.1
maf_object
str(maf_object)
maf_object@data$
maf_object@clinical.data
#Tumor_Sample_Barcode is shared in both datas

#3.1
oncoplot(maf = maf_object,
         top = 10)

ggsave("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/week7_maf/oncoplot.png")

#3.2
#The APC gene has instructions for the APC protein. The AP protein is important in tumor suppression, so a mutation in this gene could result in uncontrolled cell growth

#3.3
clinic <- maf_object@clinical.data
clinic$age_category
length(clinic$age_category)
young_mask <- clinic$age_category == "Young"
young_patient_ids <- clinic$Tumor_Sample_Barcode[young_mask]
old_patient_ids <- clinic$Tumor_Sample_Barcode[!young_mask]
length(young_patient_ids)
length(old_patient_ids)

young_maf <- subsetMaf(maf = maf_object,
                      tsb = young_patient_ids)
old_maf <- subsetMaf(maf = maf_object,
                      tsb = old_patient_ids)
#3.4
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Mutations in Young Patients", 
           m2Name = "Mutations in Old Patients")
dev.off()

ggsave("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/week7_maf/doubleoncoplot.png")

#3.5
#Most of the genes are moer mutated in older patients. In the younger patients, the genes SYNE1 and MUC16 are more mutated. This is not what I expected because I expected to see more mutations in the older patients for every gene. 

#3.6
lollipopPlot(maf_object, gene = "APC")
ggsave("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/week7_maf/lollipopplot.png")

#3.7
lollipopPlot2(m1 = young_maf,
              m2 = old_maf,
              m1_name = "Mutations in APC gene in Young Patients",
              m2_name = "Mutations in APC gene in Old Patients",
              gene = "APC")
#The gene is more commonly mutated in older patients. There are more mutations in the APC_crr region towards the end. This could be the most important part for tumor suppression. The most common mutations are nonsense mutations. There are mutations in more spots in older patients than in younger, and a lot of the mutation spots in younger patients can be found in older patients.

#3.8 
#There are 30 samples that only have both a mutation in gene b and not gene A

#3.9
#b = 7
#c = 2
#d = 35
#e = 37
#f = 42
#It does appear that A and B are not independent because their mutations appear together for most of the data.

#3.10
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")
geneB_maf <- subsetMaf(maf = maf_object,
             genes = "KRAS")

#3.11
#subsetMaf() created two objects for each gene.
geneA_maf@data$Mutation_Status
#yes, each sample only has one mutation
geneA_maf@data$Tumor_Sample_Barcode
length(geneA_maf@data$Tumor_Sample_Barcode)
#The number of samples are different. Not every sample in the original data could have had data for the TP53 gene.

#3.12
mut_bc_geneA <- geneA_maf@clinical.data$Tumor_Sample_Barcode
num_mut_geneA <- length(mut_bc_geneA) #213
mut_bc_geneB <- geneB_maf@clinical.data$Tumor_Sample_Barcode
num_mut_geneB<- length(mut_bc_geneB) #163
mut_bc_geneAB <- intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB <- length(mut_bc_geneAB)

#3.13
num_mut_geneA_only <- 135
num_mut_geneB_only <- 85

#3.14
length(maf_object@clinical.data$Tumor_Sample_Barcode)
num_neither_mutation <- 99
contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)
contig_table

#3.15
fe_results <- fisher.test(contig_table)
fe_results
#The p-value > 0.05, so the null hypothesis that the two gene mutations are independent might be false.