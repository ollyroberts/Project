library("ggplot2")
library("dplyr")
library("mixtools")

load("proline.Rda")
load("non_proline.Rda")

proline_angles <- combined_1gap_proline_16_07_19_output[,4]
#proline_angles <- proline[,3]
#non_proline_angles <- non_proline[,3]

save(proline_angles,file='proline.Rda')
save(non_proline_angles,file='non_proline.Rda')
df_proline_angles <- data.frame(proline_angles)
proline_angles <- unclass(proline_angles)

unlist(proline_angles, use.names=FALSE)

options(scipen = 999)



typeof(proline_angles)
typeof(non_proline_angles)



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





#mixmdl <- normalmixEM(proline_angles, k = 3, maxit=10000)
#summary(mixmdl)
#data.frame(bend_angle = mixmdl$x) %>%
#  
#  ggplot() +
#  geom_histogram(aes(bend_angle, ..density..), binwidth = 1, colour = "black", 
#                 fill = "white") +
#  stat_function(geom = "line", fun = plot_mix_comps,
#                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
#                colour = "red", lwd = 1.5) +
#  stat_function(geom = "line", fun = plot_mix_comps,
 #               args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
 #               colour = "blue", lwd = 1.5) +
 # stat_function(geom = "line", fun = plot_mix_comps,
#               args = list(mixmdl$mu[3], mixmdl$sigma[3], lam = mixmdl$lambda[3]),
#               colour = "green", lwd = 1.5) +
#  ylab("Density")

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
  ggtitle("two gaussian mixed model of 1Gap proline") 
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
    ggtitle("Three gaussian mixed model of 1Gap proline") 
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


