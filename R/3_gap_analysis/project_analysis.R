#install.packages("ggplot2", dependencies = TRUE)
#install.packages("mclust", dependencies = TRUE)
#install.packages("mixtools", dependencies = TRUE)
library(ggplot2)
library(mclust)
library(mixtools)

##save(proline,file="proline.Rda")
##save(non_proline,file="non_proline.Rda")

load("proline.Rda")
load("non_proline.Rda")


proline_angles <- proline[,3]
non_proline_angles <- non_proline[,3]
qplo

qplot(non_proline_angles, bin)
#mod1 = Mclust(non_proline_angles)
summary(mod1)
mod1$parameters

pro_plot <- qplot(proline_angles)
#mod2 = Mclust(proline_angles)
summary(mod2)
mod2$parameters

mod2_means <- mod2$parameters$mean ; mod2_means
pro_plot + geom_vline(xintercept = mod2_means)

#ggplot(non_proline_angles, aes(x=weight)) + 
#  geom_histogram(binwidth=1)

non_pro_angles_df <- data.frame(non_proline_angles)
geom_histogram(non_pro_angles_df, binwidth = 1)

# This is for generating mixed gaussian of second function 
geom_histogram(binwidth=1)
               
