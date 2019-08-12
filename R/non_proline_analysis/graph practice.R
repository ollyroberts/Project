library("ggplot2")
library("dplyr")
library("mixtools")
library("ggplot2")
setwd("/home/oliver/Git/Project/R/non_proline_analysis")

write(1:10,file="test.txt")
getwd()


oneGap       <- read.delim("1gap_only_03.08.2019_2angle.txt", header = FALSE, sep = "", dec = ".")
twoGap       <- read.table("2gap_only_06.08.2019.txt", header = FALSE, sep = " ", dec = ".")
threeGap     <- read.table("3_gap_only_06.08.19.txt", header = FALSE, sep = " ", dec = ".")
threeAndLess <-read.table("3gap_and_less_03.08.2019_2angle.txt", header = FALSE, sep = " ", dec = ".")
View(oneGap$V6)

breaks <- (seq(1,180,by=1))
breaksDf <- data.frame(breaks)
one_gap_angles <-data.frame(oneGap$V6)

combo <- c(oneGap$V6)



ggplot(breaks,aes(x = one_gap_angles)) +
  geom_histogram(breaks = breaks,aes(x=one_gap_angles,y=..density..),      # Histogram with density instead of count on y-axis
    col="black",fill="white") +
  geom_density(col=2) +  # Overlay with transparent density plot
  ggtitle("non proline bend angles with density function")
  


