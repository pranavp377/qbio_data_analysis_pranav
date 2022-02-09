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
  
#7 
  #I learned that females have a higher survival probability than males after a short amount of time (500 time units). Height varies, with short people showing a slightly less probabilty at the beginning, but higher probability of survival at the later times.
  
  
#Coding
  #1
  
      boxplot(c$height ~ c$gender, xlab = "Gender", ylab = "Height (cm)", main = "Height Distribution per Each Gender") #creates boxplot that compares the height distributions in females and males
      #the boxplot shows that the results are inconclusive, as there is overlap in the two boxplots. however, the average height in males is higher than that of females in the graph.
      
#2      
      install.packages("survival")
      library("survival")
      install.packages("survminer")
      library("survminer") 
      library(dplyr)
      
      clinic$days_to_death
      clinic$days_to_last_follow_up
      clinic$days_to_death <- coalesce(clinic$days_to_death, clinic$days_to_last_follow_up) #replacing NA values with corresponding values from last followup     
      clinic$days_to_death      
      
      clinic$death_event <- ifelse(clinic$vital_status == "Alive", 0, 1) #creates a column that tells the patient's current vital status
      #first variable
      surv_object <- Surv(time = clinic$days_to_death, 
                          event = clinic$death_event) #creates a survival object
      gender_fit <- surv_fit( surv_object ~ clinic$gender, data = clinic ) #creates a fit object for gender that will be used to create the survival plot
      survplot = ggsurvplot(gender_fit, 
                            pval=TRUE, 
                            ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), #plotting of the graph and creating a legend
                            legend = "right")
      p = survplot$plot + 
        theme_bw() +  # changes the appearance
        theme(axis.title = element_text(size=20), # increase font sizes
              axis.text = element_text(size=16),
              legend.title = element_text(size=14),
              legend.text = element_text(size=9))
      p
      #second variable
     clinic <- clinic[!is.na(clinic$height),]
     clinic$catheight <- ifelse(c$height < 165, "Short", "Tall") #creates categorical data for height
     surv_object <- Surv(time = clinic$days_to_death, 
                         event = clinic$death_event) #creates a survival object
     height_fit <- surv_fit( surv_object ~ clinic$catheight, data = clinic ) #creates a fit object for height that will be used to create the survival plot
     survplot = ggsurvplot(height_fit, 
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), #plotting of the graph and creating a legend
                           legend = "right") #plotting of graph and creating legend
     q = survplot$plot + 
       theme_bw() +  # changes the appearance
       theme(axis.title = element_text(size=20), # increase font sizes
             axis.text = element_text(size=16),
             legend.title = element_text(size=14),
             legend.text = element_text(size=9))
     q
      