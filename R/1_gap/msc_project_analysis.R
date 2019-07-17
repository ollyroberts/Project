library("ggplot2")
library("dplyr")
library(nortest) 

load("proline.Rda")
load("non_proline.Rda")


proline_angles <- proline[,3]
non_proline_angles <- non_proline[,3]

#skewness(non_proline_angles)

mean(non_proline_angles)
# [1] 168.9382
sd(non_proline_angles)
# [1] 14.36636
ad.test(non_proline_angles)
ks.test(non_proline_angles,pnorm(10000,168.9382,14.36636))

#normality tests 

#Two-sample Kolmogorov-Smirnov test

#data:  non_proline_angles and pnorm(10000, 168.9382, 14.36636)
#D = 1, p-value = 0.27
#alternative hypothesis: two-sided
# Ho is that both are from the same population. Ho is accepted as p_val >0.05


#
#data:  non_proline_angles
#A = 17486, p-value < 0.00000000000000022
# this indicates strongly non normal data

#ks.test(non_proline_angles)