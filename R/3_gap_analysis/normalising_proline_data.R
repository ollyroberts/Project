#  attempting to normalise the proline data by fliping values x> the maximum point. Dividing by 2 and flipping across the maximum point again
library("ggplot2")
library("dplyr")
library("mixtools")

load("proline.Rda")
load("non_proline.Rda")

#' Plot a Mixture Component
#' 
#' @param x Input data
#' @param mu Mean of component
#' @param sigma Standard deviation of component
#' @param lam Mixture weight of component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}


proline_angles <- proline[,3]

non_proline_angle_df <- data.frame(non_proline_angles)
ggplot(non_proline_angle_df, aes(x=non_proline_angles)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  coord_cartesian(xlim = c(75, 180)) +
  geom_density(alpha=.2, fill="#FF6666") +  # Overlay with transparent density plot
  ggtitle("non proline bend angles with density function")

test_vector <- c(1,2,3,4,5,6,7,8,9,10)
test_df <- data.frame(test_vector)

for (i in test_vector){
  if (test_vector[i] > 5){
    difference<- test_vector[i] - 5
    test_vector[i] = test_vector[i] - (difference*2)
  }
}
greater_than_5 <- test_df[test_df>5]