library("ggplot2")
#library("dplyr")
library("mixtools")

library(fitdistrplus)
setwd("/home/oliver//Git/Project/R/3_gap_analysis/")
load("3_gap_proline.Rda")
load("non_proline.Rda")




options(scipen = 999)

load("~/Git/Project/R/3_gap_analysis/3_gap_proline.Rda")
load("~/Git/Project/R/3_gap_analysis/non_proline.Rda")

proline_angles <- proline[,3]
non_proline_angles <- non_proline[,3]
non_proline_angleDF <- data.frame(non_proline_angles)
proline_angleDF <- data.frame(proline_angles)


# Fitting distributions with MAss or fitdistplus 
fit <- fitdistr(non_proline_angles, densfun="normal")  # we assume my_data ~ Normal(?,?)
hist(non_proline_angles, pch=20, breaks=.5, prob=TRUE, main="")
curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
print(fit$loglik)
print(fit$mean)
summary(fit$estimate$Mean)
attributes(fit$estimate)

ggplot(non_proline_angleDF, aes(X = non_proline_angles)) +
  geom_histogram(aes(non_proline_angles, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean =171.94 , sd = 4.2 ), colour = "Red", lwd = 1.5) +
#stat_function(aes(x =non_proline_angles, y =..density..),fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) + ylab("") +
#  scale_y_continuous(breaks = NULL)  
#  stat_function(geom = "line", fun = dnorm, fun
#                args = list(168.93817179, 14.36630247, lam = 1),
#                colour = "red", lwd = 1.5) + 
  coord_cartesian(xlim = c(75, 180)) +
ylab("Density") +
xlab(" Bend Angle")

# fitting a single curve to 3_gap proline data
fit <- fitdistr(proline_angles, densfun="normal")  # we assume my_data ~ Normal(?,?)
summary(fit)

ggplot(proline_angleDF, aes(X = proline_angles)) +
  geom_histogram(aes(proline_angles, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(fun = dnorm, args = list(mean =168.94 , sd = 14.3663 ), colour = "Red", lwd = 1.5) +
  #stat_function(aes(x =non_proline_angles, y =..density..),fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) + ylab("") +
  #  scale_y_continuous(breaks = NULL)  
  #  stat_function(geom = "line", fun = dnorm, fun
  #                args = list(168.93817179, 14.36630247, lam = 1),
  #                colour = "red", lwd = 1.5) + 
  coord_cartesian(xlim = c(75, 180)) +
  ylab("Density") +
  xlab(" Bend Angle")

#' Plot a Mixture Component
#' 
#' @param x Input data
#' @param mu Mean of component
#' @param sigma Standard deviation of component
#' @param lam Mixture weight of component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)


mixmdl2 <- (normalmixEM(proline_angles, k = 2, maxit=10000))
mixmdl3 <- (normalmixEM(proline_angles, k = 3, maxit=10000))
#mixmdl4 <- (normalmixEM(proline_angles, k = 4, maxit=10000))
#mixmdl5 <- (normalmixEM(proline_angles, k = 5, maxit=10000))
#mixmdl6 <- (normalmixEM(proline_angles, k = 6, maxit=10000))
#mixmdl7 <- (normalmixEM(proline_angles, k = 7, maxit=10000))





mixmdl <- normalmixEM(non_proline_angles,lambda = 1, maxit=10000)
summary(mixmdl)
#data.frame(bend_angle = mixmdl$x) %>%
#  
  ggplot() +
  geom_histogram(aes(bend_angle, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
 # stat_function(geom = "line", fun = plot_mix_comps,
#               args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
#               colour = "green", lwd = 1.5) +
  coord_cartesian(xlim = c(75, 180)) +
  ggtitle("Two gaussian mixed model for 'Normal Helix' data ") 
  ylab("Density")

summary(mixmdl2)
data.frame(bend_angle = mixmdl2$x) %>%
  
  ggplot() +
  geom_histogram(aes(bend_angle, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl2$mu[1], mixmdl2$sigma[1], lam = mixmdl2$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl2$mu[2], mixmdl2$sigma[2], lam = mixmdl2$lambda[2]),
                colour = "blue", lwd = 1.5) +
  coord_cartesian(xlim = c(75, 180)) +
  ggtitle("two gaussian mixed model") 
  ylab("Density")

summary(mixmdl3)
  data.frame(bend_angle = mixmdl3$x) %>%
    
    ggplot() +
    geom_histogram(aes(bend_angle, ..density..), binwidth = 1, colour = "black", 
                   fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl3$mu[1], mixmdl3$sigma[1], lam = mixmdl3$lambda[1]),
                  colour = "red", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl3$mu[2], mixmdl3$sigma[2], lam = mixmdl3$lambda[2]),
                  colour = "blue", lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl3$mu[3], mixmdl3$sigma[3], lam = mixmdl3$lambda[2]),
                  colour = "green", lwd = 1.5) +
    coord_cartesian(xlim = c(75, 180)) +
    ggtitle("Three gaussian mixed model") 
  ylab("Density")

proline_angle_df <- data.frame(proline_angles)
data.frame(bend_angle = mixmdl3$x) %>%
ggplot(proline_angle_df, aes(x=proline_angles)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  coord_cartesian(xlim = c(75, 180)) +

  geom_density(alpha=.2, fill="#FF6666") +  # Overlay with transparent density plot
ggtitle("proline bend angles with density function")


non_proline_angle_df <- data.frame(non_proline_angles)
ggplot(non_proline_angle_df, aes(x=non_proline_angles)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  coord_cartesian(xlim = c(75, 180)) +
  geom_density(alpha=.2, fill="#FF6666") +  # Overlay with transparent density plot
  ggtitle("non proline bend angles with density function")


