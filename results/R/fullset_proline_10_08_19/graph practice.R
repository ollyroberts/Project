library("ggplot2")
library("dplyr")
library("mixtools")
library("ggplot2")
library(readr)



setwd("/home/oliver/Git/Project/R/fullset_proline_10_08_19")
#write(1:10,file="test.txt")
getwd()
write.csv(combo1_2_3, file = "one_two_three_gap_prolines.csv")



breaks <- (seq(1,180,by=1))
breaksDf <- data.frame(breaks)

one_gap_angles <-oneGap$V6


no_gap_H <- data.frame(read_table2("no_gap_H_12_08_19.txt", 
                                 col_names = FALSE))
no_gap_h <- data.frame(read_table2("no_gap_h_12_08_19.txt", 
                                col_names = FALSE))
one_gap_not_horH <- data.frame(read_table2("one_gap_not_horH_12_08_19.txt", 
                                col_names = FALSE))
two_gap<- data.frame(read_table2("two_gap_12_08_19.txt", 
                                col_names = FALSE))
three_gap <- data.frame(read_table2("three_gap_12_08_19.txt", 
                                col_names = FALSE))

no_gap_H$gap  <- "0"
no_gap_h$gap <-  "0"
one_gap_not_horH$gap <- "1"
two_gap$gap <- "2"
three_gap$gap <- "3"

no_gap_HDF <- data.frame(c(no_gap_H["X1"],no_gap_H["X6"],no_gap_H["X10"],no_gap_H["gap"]))
no_gap_hDF <- data.frame(c(no_gap_h["X1"],no_gap_h["X6"],no_gap_h["X10"],no_gap_h["gap"]))
one_gap_not_horHDF <- data.frame(c(one_gap_not_horH["X1"],one_gap_not_horH["X6"],one_gap_not_horH["X10"],one_gap_not_horH["gap"]))
two_gapDF <- data.frame(c(two_gap["X1"],two_gap["X6"],two_gap["X10"],two_gap["gap"]))
three_gapDF <- data.frame(c(three_gap["X1"],three_gap["X6"],three_gap["X10"],three_gap["gap"]))

combined <- rbind(no_gap_HDF,no_gap_hDF,one_gap_not_horHDF,two_gapDF,three_gapDF)
names(combined)[names(combined) == "X6"] <- "angle"
names(combined)[names(combined) == "X10"] <- "secstr"
names(combined)[names(combined) == "X1"] <- "pdb"

only_no_gap <- rbind(no_gap_H,no_gap_h)
names(only_no_gap)[names(only_no_gap) == "X6"] <- "angle"
names(only_no_gap)[names(only_no_gap) == "X10"] <- "secstr"
names(only_no_gap)[names(only_no_gap) == "X1"] <- "pdb"

# To remove a row using criteria
combined_by_secstr <- combined[!(combined$secstr=="NA"),]
#combined_by_secstr <- combined_by_secstr[!(combined_by_secstr$secstr=="S"),]
combined_no_na <-combined_by_secstr[complete.cases(combined_by_secstr), ]
table(combined$gap)


only_no_gap_secstr <- 

combined_by_secstr <- data[(combined$angle =="H"),]
sum(combined_no_na$angle > 164.48)
# Remove NA groups from no_gaph and no_gap_H


only_no_gap_secstr <- subset(only_no_gap, secstr=="H" | secstr=="h")

only_no_gap_secstr <- subset(only_no_gap, angle<=164.48)
only_no_gap_secstr_H_164 <- subset(only_no_gap, secstr=="H" & angle>=164.48)
only_no_gap_secstr_h_164 <- subset(only_no_gap, secstr=="h" & angle>=164.48)

sum(only_no_gap_secstr_H_164$angle > 164.48)
sum(only_no_gap_secstr_h_164$angle  > 164.48)

table(only_no_gap_secstr$secstr)

sum

 #one_gap_not_horH$type <- data.frame(one_gap_not_horH)
#two_gap$type  <- data.frame(two_gap)
#three_gap$type <- data.frame(three_gap)

#comb1 <- merge(no_gap_HDF,no_gap_hDF)

# import the files and asign the angle vector (in this case column 6) to a dataframe
#no_gap_anglesDf <-
#one_gap_anglesDf   <- data.frame(X1gap_only_03_08_2019_2angle$X6)
#two_gap_anglesDf   <- data.frame(X2gap_only_06_08_2019$X6)
#three_gap_anglesDf <- data.frame(X3_gap_only_06_08_19$X6)

#one_gap_anglesDf$Type   <- "1"
#two_gap_anglesDf$Type   <- "2"
#three_gap_anglesDf$Type <- "3"

#names(one_gap_anglesDf)[names(one_gap_anglesDf) == "X1gap_only_03_08_2019_2angle.X6"] <-"angles"
#names(two_gap_anglesDf)[names(two_gap_anglesDf) == "X2gap_only_06_08_2019.X6"] <-"angles"
#names(three_gap_anglesDf)[names(three_gap_anglesDf) == "X3_gap_only_06_08_19.X6"] <-"angles"


#combo1_2 <-rbind(one_gap_anglesDf,two_gap_anglesDf)
#combo1_2_3 <-rbind(combo1_2,three_gap_anglesDf)





oneGapNumeric <- numeric(oneGap$V6)

ggplot(oneGap,aes(x = oneGap$V6)) +
  geom_histogram(alpha =0.2,breaks = breaks,aes(x=oneGap$V6))+
  geom_density(aes(x=oneGap$V6),col="red")
 
# using the type ggplot(data, aes(x = rating, fill = type)) to seperate out data 
# to use this we must combine ooneGap,twoGap,ThreeGap
  
#ggplot(oneGap,aes(x = oneGap$V6)) +
#  geom_density() +
#  geom_density(aes(x=twoGap$V6,col="red")) +
#  geom_histogram(alpha = 0.1, aes(y = ..density..))


#ggplot(combined, aes(x = angle, fill = gap)) +
  ##geom_histogram(binwidth = .5, alpha =.5, position = "identity") +
  #geom_density(alpha = .3)+
  #ggtitle("Density of Proline bend angle  in alpha helices seperated by gap length") 
  #ylab("Density") +
  #xlab("Bend angle of helices")


# Add weighting for combined 
combined <- combined %>%
  group_by(secstr) %>%
  mutate(n = n()/nrow(combined))

ggplot(combined, aes(x=angle,fill=secstr)) +
  geom_density(aes(weight = n), col=NA, alpha = 0.35)

ggplot(combined, aes(x = secstr, y = angle, fill=secstr))+
  geom_violin( col = NA)

# Chosen  
ggplot(combined, aes(x = angle, fill = gap)) +
  geom_histogram(binwidth = 1, alpha =.5, position = "identity") +
  #geom_density(alpha = .3)+
  coord_cartesian(xlim = c(75, 180)) +
  ggtitle("histogram of 'gap' containing Proline helices in alpha helices grouped by gap length") 
  ylab("Density") +
  xlab("Bend angle of helices") +
  
  
# This is an attempt to add density curve to diagram
ggplot(combined, aes(x = angle, fill = gap)) +
  geom_density(alpha = .3)+
  ggtitle("histogram of 'gap' containing Proline helices in alpha helices seperated by gap length") +
  ylab("Density") +
  xlab("Bend angle of helices")


ggplot(combined, aes(x = angle, fill = secstr)) +
  #geom_histogram(binwidth = .5, alpha =.5, position = "identity") +
  geom_density(alpha = .3)+
  ggtitle("Density of all Proline bend angles in alpha helices seperated by Secondary structure, this includes") 
  ylab("Density") +
  xlab("Bend angle of helices")

# Chosen 
ggplot(combined_by_secstr, aes(x = angle, fill = secstr)) +
  geom_histogram(binwidth = 1, alpha =.5, position = "identity") +
  #geom_density(alpha = .3)+
  #ggtitle("Histogram of all 1-3 Proline bend angles in alpha helices seperated by Secondary structure, this includes") 
  ylab("Density") +
  xlab("Bend angle of helices")  

  
ggplot(only_no_gap, aes(x = angle, fill = secstr)) +
  #geom_histogram(binwidth = .5, alpha =.5, position = "identity") +
  geom_density(alpha = .3)+
  #ggtitle("Density of Proline bend angles in alpha helices with no gap, seperated by secondary structure") 
  ylab("Density") +
  scale_x_continuous(name="Bend angle of helices")

only_no_gap
only_no_gap_secstr <- subset(only_no_gap, secstr=="H" | secstr=="h")
# chose density of no gap  
ggplot(only_no_gap_secstr, aes(x = angle, fill = secstr)) +
  geom_histogram(binwidth = 0.5, alpha =0.3, position = "identity") +
  #geom_density(alpha = .3)+
  #ggtitle("Histogram of Proline bend angles in alpha helices with no gap, seperated by secondary structure") 
  ylab("Count") +
  xlab("Bend angle of 'no gap' helices")

# Maybe a violin plot or weighted geom_density is the solution https://www.youtube.com/watch?v=-jO2wuu7GUQ
# first we must add weights 
only_no_gap <- only_no_gap %>%
  group_by(secstr) %>%
  mutate(n = n()/nrow(only_no_gap))

# modify only_no_gap to remove Na
  ggplot(only_no_gap, aes(x=angle,fill=secstr)) +
    geom_density(aes(weight = n), col=NA, alpha = 0.1)
  
ggplot(only_no_gap, aes(x = secstr, y = angle, fill=secstr))+
  geom_violin(aes(weight =(0.91458634,0.07107489)), col = NA)

