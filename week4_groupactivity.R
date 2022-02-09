# this will install packages (if necessary) and load them
if(!require(BiocManager)) install.packages("BiocManager")

# the double colon syntax calls a function from a specific package
# this avoids loading the entire package
# in this case, we need to download TCGAbiolinks from Bioconductor using BiocManager
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

# this just loads a package
library(TCGAbiolinks)
clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
# Only use this line ONCE! Comment out after you have downloaded the data. 
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

clinic <- read.csv("clinic.csv")

#Written Activity
#1
#  A categorical variable is a variable that can be put into distinct, finite categories. For example, hair color is a variable that can be divided into categories such as brown, black, blonde, etc.
# A discrete variable is a variable that is collected through counting. For example, the number of people in a classroom is a discrete variable.
# A continuous variable is a variable that is collected through measuring. For example, the height of people in a classroom is a continuous variable because it is measured.

#2
  clinic$gender
  is.na(clinic$gender)
# I have chosen the gender variable. There are 0 NA values in this column.
  
#3 
# Gender can be measured/collected through asking the participant/patient. Gender is a categorical variable, as it has a number of categories that the data can be organized into.
  
#4
  # https://www.frontiersin.org/articles/10.3389/fonc.2020.607909/full#:~:text=A%20higher%20incidence%20of%20colorectal,in%20CRC%20rates%20and%20survival.
    #The study is titled "Sexual Dimorphism in Colon Cancer."  The results suggest that there is a higher incidence of CRC in males than in females. Young women also might have a better survival outcome than men in CRC cases. It is theorized that the hormone estrogen could play a role in these results
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4419057/
    #The study is titled "Sex- and gender-specific disparities in colorectal cancer risk." This study explores the statistic that females over 65 have a lower survival rate than men. The study highlights many flaws within the treatment of CRC, in which gender is often not considered when treatment is administered to patient. There is often a difference in tumour location in men and women, so sex and gender are important factors to consider. 

#5 
  clinic$height
  is.na(clinic$height) #coutning NA values
  c = clinic[!is.na(clinic$height),] #removing NA values
  c$height #makings sure NA values are remoevd
  #Height is a continuous variable as it is measured. This variable can be collected through the measuring of a person's height.
  
#6 
  #The gender data collected in this dataset was stored either as "Male" or "Female."
  # 1) Hypothesis 1: Patients that are categorized as male are likely to have a greater average height than patients categorized as female.
  # 2) Hypothesis 2: The survival in colorectal cancer of patients categorized as female is likely to be higher than that of patients categorized as male.
  # 3) Hypothesis 3: Patients with a lower height are likely to have a higher survival than patients with a higher height.
  
  
#Coding
  #1
  
      boxplot(c$height ~ c$gender, xlab = "Gender", ylab = "Height (cm)", main = "Height Distribution per Each Gender") #creates boxplot that compares the height distributions in females and males
      #the boxplot shows that the results are inconclusive, as there is overlap in the two boxplots. however, the average height in males is higher than that of females in the graph.
  
  