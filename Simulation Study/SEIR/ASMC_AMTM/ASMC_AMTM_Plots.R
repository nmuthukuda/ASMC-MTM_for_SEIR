setwd("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//ASMC_AMTM")
# Loading required packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(coda)
library(gridExtra)
library(purrr)

resampled_particles1 <- read.csv("ASMC_AMTM_try_rsp.csv")

mcmc_object1 <- mcmc(resampled_particles1)

# Summary 
summary(mcmc_object1)

resampled_particles1 <- resampled_particles1 %>%
  dplyr::rename(Beta = V1, Delta = V2, Gamma = V3, Phi = V4, S0 = V5, E0 = V6, I0 = V7)


# Reshape the data to long format
long_data <- resampled_particles1 %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

# Plot density curves with facets
p1 <- ggplot(long_data, aes(x = Value, fill = Parameter)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Parameter, scales = "free") +
  labs(title = "",
       x = "Value",
       y = "Density") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14), 
        strip.text = element_text(size = 14))  # Increase facet label size

# Save the plot with high resolution
ggsave("ASMC_AMTM_DEN_500_New.jpeg", plot = p1, width = 12, height = 8, dpi = 300)

# plot for ASMC_MTM_SEIR

par(mfrow = c(2,1))
ess7 <- read.csv("ASMC_AMTM_Finale_ESS.csv")
alpha7 <- read.csv("ASMC_AMTM_Finale_alpha.csv")

ess7 <- ess7 %>%
  mutate(Iteration = 1:56)

ess7 <- ess7 %>%
  dplyr::select(Iteration, everything())

alpha7 <- alpha7 %>%
  mutate(Iteration = 1:57)

alpha7 <- alpha7 %>%
  dplyr::select(Iteration, everything())

# Plotting the annealing parameter
p1 <- ggplot(alpha7, aes(x = Iteration, y = alpha7[, 2])) +
  geom_point(color = "black", fill = "black", shape = 21) +
  labs(x = 'Iteration', y = expression(alpha)) +
  ggtitle("") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16))

# Plotting the ESS values
p2 <- ggplot(ess7, aes(x = Iteration, y = ess7[, 2])) +
  geom_point(color = "black", fill = "black", shape = 16) +
  geom_line(linetype = "dashed", color = "grey") +
  labs(x = 'Iteration', y = 'ESS') +
  ggtitle("") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16))

# Arrange the plots side by side
cplot1 <- grid.arrange(p1, p2, nrow = 1)

# Save the plot with high resolution
ggsave("ASMC_AMTM_SEIR_ESS_Alpha_New.jpeg", plot = cplot1, width = 12, height = 8, dpi = 300)