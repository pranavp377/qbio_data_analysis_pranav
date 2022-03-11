library(maftools)
library(TCGAbiolinks)
library(ggplot2)

clinic <- data.table::fread("/Users/mailp/Desktop/USC/QBIODataAnalysis/qbio_data_analysis_pranav/analysis_data/coad_clinical_data.csv",
                            data.table = F)
clinic$height

#splitting clinic dataset into male and female datasets based on the gender column. 
clinic$gender
length(clinic$gender) #524
female_clinic <- clinic[ clinic$gender == 'FEMALE', ]
female_clinic$gender
length(female_clinic$gender) #244
male_clinic <- clinic[clinic$gender == 'MALE', ]
length(male_clinic$gender) #280
male_clinic$gender

#remove NA values from height column in female clinic dataset
female_clinic$height
is.na(female_clinic$height)
female_mask <- is.na(female_clinic$height)
sum(female_mask) #125 NA Values
female_clinic_new <- female_clinic[female_mask == 'FALSE', ]
female_clinic_new$height
length(female_clinic_new$height) #119

#remove NA values from height column in male clinic dataset
male_clinic$height
is.na(male_clinic$height)
male_mask <- is.na(male_clinic$height)
sum(male_mask) #147 NA Values
male_clinic_new <- male_clinic[male_mask == 'FALSE', ]
length(male_clinic_new$height) #133

#cleaned datasets
female_clinic_new
male_clinic_new




library(survival)
library(survminer)

#survival analysis for height in females
female_clinic_new$height_category <- ifelse(female_clinic_new$height >= 162.5, "Tall", "Short")
surv_object_female <- Surv(time = female_clinic_new$days_to_death,
                           event = female_clinic_new$death_event)
height_fit_female <- surv_fit(surv_object_female ~ female_clinic_new$height_category, data = female_clinic_new)

survplot = ggsurvplot(height_fit_female,
                      pval = TRUE,
                      ggtheme = theme(plot.margin = unit(c(25,25,25,25), "cm")),
                      legend = "right")

q = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=10), 
        axis.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
q

#survival analysis for males
male_clinic_new$height_category <- ifelse(male_clinic_new$height >= 175.3, "Tall", "Short")
surv_object_male <- Surv(time = male_clinic_new$days_to_death,
                           event = male_clinic_new$death_event)
height_fit_male <- surv_fit(surv_object_male ~ male_clinic_new$height_category, data = male_clinic_new)

survplot = ggsurvplot(height_fit_male,
                      pval = TRUE,
                      ggtheme = theme(plot.margin = unit(c(25,25,25,25), "cm")),
                      legend = "right")


p = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=10), 
        axis.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
p



#mutation analysis in females
tall_female_mask <- female_clinic_new$height_category == "Tall"
short_female_mask <- female_clinic_new$height_category == "Short"
colnames(female_clinic_new)[ colnames(female_clinic_new) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
tall_female_IDs <- female_clinic_new$Tumor_Sample_Barcode[tall_female_mask]
short_female_IDs <- female_clinic_new$Tumor_Sample_Barcode[short_female_mask]
length(tall_female_IDs)
length(short_female_IDs)

mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object_female <- read.maf(maf = mutation_query, 
                       clinicalData = female_clinic_new, 
                       isTCGA = TRUE)

tall_maf_female <- subsetMaf(maf = maf_object_female,
                              tsb = tall_female_IDs)
short_maf_female <- subsetMaf(maf = maf_object_female,
                              tsb = short_female_IDs,)

coOncoplot(m1 = tall_maf_female, 
           m2 = short_maf_female, 
           m1Name = "Tall Females (>= 162.5 cm)", 
           m2Name = "Short Females (< 162.5 cm)")
lollipopPlot2(m1 = tall_maf_female, 
              m2 = short_maf_female, 
              m1_name = "Tall Female Patients (>= 162.5 cm)",
              m2_name = "Short Female Patients (< 162.5 cm)",
              gene = "IGF1")

#mutation analysis in males
tall_male_mask <- male_clinic_new$height_category == "Tall"
short_male_mask <- male_clinic_new$height_category == "Short"
colnames(male_clinic_new)[ colnames(male_clinic_new) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
tall_male_IDs <- male_clinic_new$Tumor_Sample_Barcode[tall_male_mask]
short_male_IDs <- male_clinic_new$Tumor_Sample_Barcode[short_male_mask]
length(tall_male_IDs)
length(short_male_IDs)


maf_object_male <- read.maf(maf = mutation_query, 
                              clinicalData = male_clinic_new, 
                              isTCGA = TRUE)

tall_maf_male <- subsetMaf(maf = maf_object_male,
                             tsb = tall_male_IDs)
short_maf_male <- subsetMaf(maf = maf_object_male,
                              tsb = short_male_IDs,)

coOncoplot(m1 = tall_maf_male, 
           m2 = short_maf_male, 
           m1Name = "Tall Males (>= 175.3 cm)", 
           m2Name = "Short Males (< 175.3 cm)")
lollipopPlot2(m1 = tall_maf_male, 
              m2 = short_maf_male, 
              m1_name = "Tall Male Patients (>= 175.3 cm)",
              m2_name = "Short Male Patients (< 175.3 cm)",
              gene = "IGF1")

#RNASeq
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts") 
GDCdownload(query)
sum_exp <- GDCprepare(query)
counts = assays(sum_exp)$"HTSeq - Counts"
assays(sum_exp)$"HTSeq - Counts"
gene_id_mask = (rowData(sum_exp)$external_gene_name == "IGF1")
sum(gene_id_mask)
ensembl1_gene = rowData(sum_exp)$ensembl_gene_id[gene_id_mask] 
ensembl1_gene

bool_gender_na = is.na(colData(sum_exp)$gender)
num_na = sum(bool_gender_na)
num_na
gender_no_NAs = colData(sum_exp)$gender[!bool_gender_na]
length(gender_no_NAs)


library(ggplot2)

identical( rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")  ) #TRUE
gene_counts = assays(sum_exp)$"HTSeq - Counts"["ENSG00000017427",!bool_gender_na]

length(gene_counts) == length(gender_no_NAs) #TRUE  
#gene_counts was created by using a boolean mask to isolate the columns without NA values, that corresponded with the row of the gene name A

boxplot(gene_counts ~ gender_no_NAs,
        xlab = "Sex of patient",
        ylab = "Gene Counts")

