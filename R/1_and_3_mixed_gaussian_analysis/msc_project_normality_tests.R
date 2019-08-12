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



#normality tests 

ad.test(non_proline_angles)
# p-value < 0.00000000000000022
# not_normal distribution 

ks.test(non_proline_angles,pnorm(10000,168.9382,14.36636))
#Two-sample Kolmogorov-Smirnov test
#data:  non_proline_angles and pnorm(10000, 168.9382, 14.36636)
#D = 1, p-value = 0.27
#alternative hypothesis: two-sided
# Ho is that both are from the same population. Ho is accepted as p_val >0.05

#Shapiro Wilks test 
non_proline_p_values <- numeric()
for (i in 1:5000){
  temp_samp <- sample(5000,non_proline_angles)
  test_result <- shapiro.test(temp_samp)
  non_proline_p_values <- c(non_proline_p_values,test_result$p.value)
}
non_proline_p_values[0.05<non_proline_p_values]
# no pvalues >0.05 = 0
# Ho is that the data is normally distributed. Ho is rejected.

# comparison tests between proline 3rd mixed gaussian model(of 3)
# and non_proline data

# proline_3rd mixed gaussian 
# mean =170.130773 sd = 4.571051
proline_3rd_mixed_sample <- rnorm(10000,mean = 170.130773,sd = 4.571051)


# non_proline bend angles 
# mean = 168.9382 sd = 14.3663

# this means we shall use unquel varience test

t.test(proline_3rd_mixed_sample,non_proline_angles)
# Ho is that the difference in means is = to zero. 
# p-value <0.00000000000000022
# we reject Ho
