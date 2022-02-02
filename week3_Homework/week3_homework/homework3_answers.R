#exercise 1.1
attenu
attenu$station
which(is.na(attenu$station)) #every row that makes is.na() TRUE is displayed with the which() function
attenu_cleaned <- na.omit(attenu) #na.omit displays all rows without an NA value
head(attenu_cleaned)
dim(attenu_cleaned)

#exercise 1.2
Theoph
Theoph_2 <- Theoph
str(Theoph) #Dose information in column titled Dose
median(Theoph_2$Dose)
Theoph_2$Dose_Class <- ifelse(Theoph_2$Dose >= 4.53, "high", "low") 
head(Theoph_2) 

#exercise 1.3
starbucks <- read.csv("starbucks.csv")
is.na(starbucks)
is_row_empty <- vector(mode = "logical", nrow(starbucks))
nrow(starbucks)
length(is_row_empty)
starbucks_cleaned <- starbucks[(rowSums(is.na(starbucks)) < 6),]
starbucks_cleaned
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carb (g)", ylab = "Calories", col = ifelse(starbucks_cleaned$Calories == max(starbucks_cleaned$Calories), "red", "black")) #more carbs seems to correlate with more calories
max(starbucks_cleaned$Calories)
starbucks_cleaned[starbucks_cleaned$Calories == 430,]
starbucks_cleaned$is_highest_fat <- ifelse(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat), "TRUE", "FALSE")
starbucks_cleaned

#exercise 1.4
Batting <- read.csv("Batting.csv")
Batting[Batting$HR >= 3,]
plot(Batting$yearID, Batting$HR, xlab = "Year", ylab = "Number of Homeruns")
LAAngels <- Batting[Batting$teamID == "LAA", ]
LAAngels
plot(LAAngels$yearID, LAAngels$HR, xlab = "Year", ylab = "Number of Homeruns")
NewBatting <- Batting[Batting$teamID == "ATL" | Batting$teamID == "PIT",]
NewBatting$color <- ifelse(NewBatting$teamID == "ATL", "red", "blue")
plot(NewBatting$yearID, NewBatting$HR, col = NewBatting$color, xlab = "Year", ylab = "Number of Homeruns")

#exercise 1.5
easy_plot <- function(x,y,color_data) {
  m <- median(color_data)
  levels <- ifelse(color_data < m, "low", "high")
  levels = factor(levels)
  return(m)
  plot(x,y,pch=20,col = levels)
  cor.test(x,y)

}

#exercise 2.1
#The iris data set describes the sepal length, sepal width, petal length, and petal width of three different species of flower. For each datapoint, 5 features are present (sepal length, sepal width, petal length, petal width, species name)
#exercise 2.2
# The Sepal.Length, Sepal.Width, Petal.Length, and Petal.Width data points are continuous. The species data points are factors. 
class(iris$Sepal.Length) #numeric
class(iris$Sepal.Width) #numeric
class(iris$Petal.Length) #numeric
class(iris$Petal.Width) #numeric
class(iris$Species) #factor
#exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)
#The petal data points appear to both be skewed to the right a bit. The sepal datasets, especially the width, follow a more symmetric structure.
#Exercise 2.4
iris_copy <- iris
meansepalw <- mean(iris$Sepal.Width)
meansepalw
iris$Sepal.Size <- ifelse(iris$Sepal.Width <= 3.057333, "narrow", "wide")
boxplot(iris$Sepal.Width ~ iris$Sepal.Size, xlab = "Sepal Size", ylab = "Sepal Width")
#Excercise 2.5
#The species represented by the color black looks to be the most unique. The species represented by green and orange look the most similar.
iris_copy <- subset(iris_copy, select = -Species)
pairs(iris_copy)
install.packages("GGally")
library("GGally")
ggpairs(iris_copy)

#Exercise 3.1
install.packages("TCGAbiolinks")
library("TCGAbiolinks")
