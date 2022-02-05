if(!require(BiocManager)) install.packages("BiocManager")
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml") 
GDCdownload(clin_query) #data downloaded
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
#Excercise 1.1
head(clinic)
str(clinic)

#1.2
colnames(clinic)
clinic$informed_consent_verified


clinic$race_list = as.character(clinic$race_list)

#Exercise 2.1
clinic$race_list
unique(clinic$race_list)
clinic$age_at_initial_pathologic_diagnosis
#2.2
boxplot(clinic$race_list ~ clinic$age_at_initial_pathologic_diagnosis, xlab = "Age at initial diagnosis", ylab = "Race", las = 2, cex.axis = 0.5)
par(mar=c(10,1,1,1))

#2.3
sum(clinic$race_list == "") #sums the number of columns with an empty value
replace(clinic$race_list, clinic$race_list == "", "No data")
library(tidyr)

clinic$race_list <- clinic$race_list %>% replace_na('No data')

#2.4
clinic$age_at_initial_pathologic_diagnosis
summary(clinic$age_at_initial_pathologic_diagnosis)

#2.5
# young <- 0
#old <- 0
#length(clinic$age_at_initial_pathologic_diagnosis)
#for (i in 1:524) {
  
#  if (clinic$age_at_initial_pathologic_diagnosis[i] < 50) {
#    young <- young + 1
    
#  } else (clinic$age_at_initial_pathologic_diagnosis[i] >= 50) {
#    old <- old + 1
#  }
#  return(old)
#  return(young)
#} #this code keeps resulting in an error "unexpected }" 

sum(clinic$age_at_initial_pathologic_diagnosis < 50)
sum(clinic$age_at_initial_pathologic_diagnosis >= 50)
length(clinic$age_at_initial_pathologic_diagnosis)

#2.6
clinic$bcr_patient_barcode

clinic$age_category <- ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "Young", "Old" )
young.mask <- clinic$age_category == "Young"
clinic[clinic$bcr_patient_barcode,]

#2.8
clinic[1,1] #1st row, 1st column
clinic[1,] #first row, every column
clinic[2:5,] #second, third, fouth, fifth row and every column
clinic[,3] #every row, third column







#3
install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)
library(dplyr)

#3.1
clinic$days_to_death
clinic$days_to_last_followup
clinic$days_to_death = coalesce(clinic$days_to_death, clinic$days_to_last_followup)
clinic$days_to_death
#3.2
clinic$death_event <- ifelse(clinic$vital_status == "Alive", 2, 1 )

clinic$death_event
clinic$vital_status
 
#3.3
surv_object <- Surv(time = clinic$days_to_death, event = clinic$death_event)
race_fit <- surv_fit(surv_object ~ clinic$race_list, data = clinic)
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=5))
p
ggsave("../week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)
