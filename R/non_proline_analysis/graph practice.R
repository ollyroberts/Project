library("ggplot2")
library("dplyr")
library("mixtools")
library("ggplot2")
setwd("/home/oliver/Git/Project/R/non_proline_analysis")

write(1:10,file="test.txt")
getwd()
write.csv(combo1_2_3, file = "one_two_three_gap_prolines.csv")



breaks <- (seq(1,180,by=1))
breaksDf <- data.frame(breaks)



one_gap_angles <-oneGap$V6


# import the files and asign the angle vector (in this case column 6) to a dataframe
one_gap_anglesDf   <- data.frame(X1gap_only_03_08_2019_2angle$X6)
two_gap_anglesDf   <- data.frame(X2gap_only_06_08_2019$X6)
three_gap_anglesDf <- data.frame(X3_gap_only_06_08_19$X6)

one_gap_anglesDf$Type   <- "1"
two_gap_anglesDf$Type   <- "2"
three_gap_anglesDf$Type <- "3"

names(one_gap_anglesDf)[names(one_gap_anglesDf) == "X1gap_only_03_08_2019_2angle.X6"] <-"angles"
names(two_gap_anglesDf)[names(two_gap_anglesDf) == "X2gap_only_06_08_2019.X6"] <-"angles"
names(three_gap_anglesDf)[names(three_gap_anglesDf) == "X3_gap_only_06_08_19.X6"] <-"angles"


combo1_2 <-rbind(one_gap_anglesDf,two_gap_anglesDf)
combo1_2_3 <-rbind(combo1_2,three_gap_anglesDf)





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


ggplot(combo1_2_3, aes(x = angles, fill = Type)) +
  #geom_histogram(binwidth = .5, alpha =.5, position = "identity") +
  geom_density(alpha = .3)+
  ggtitle("Density of Proline bend anglein alpha helices seperated by gap length") 
  ylab("Density") +
  xlab("Bend angle of helices")






