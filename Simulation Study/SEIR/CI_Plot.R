setwd("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR")
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the CSV files
chain1 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC//AMCMC_Finale_chain1.csv")
chain2 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC//AMCMC_Finale_chain2.csv")
chain3 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC//AMCMC_Finale_chain3.csv")
chain4 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC//AMCMC_Finale_chain4.csv")
chain5 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC//AMCMC_Finale_chain5.csv")

chain11 <- chain1[25001: 59000, ]
chain22 <- chain2[25001: 59000, ]
chain33 <- chain3[25001: 59000, ]
chain44 <- chain4[25001: 59000, ]
chain55 <- chain5[25001: 59000, ]


library(coda)
mcmc_chain1 <- mcmc(chain1)
mcmc_chain2 <- mcmc(chain2)
mcmc_chain3 <- mcmc(chain3)
mcmc_chain4 <- mcmc(chain4)
mcmc_chain5 <- mcmc(chain5)

mcmc_list <- mcmc.list(mcmc_chain1, mcmc_chain2, mcmc_chain3, mcmc_chain4, mcmc_chain5)

summary(mcmc_list)

resampled_particles1 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//ASMC_AMTM//ASMC_AMTM_try_rsp.csv")

mcmc_object1 <- mcmc(resampled_particles1)

# Summary 
summary(mcmc_object1)

resampled_particles2 <- read.csv("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AASMC//AASMC_Finale_1000_0.99_rsp.csv")

mcmc_object2 <- mcmc(resampled_particles2)

# Summary 
summary(mcmc_object2)


# Calculate the 95% HPD intervals for the combined MCMC chains
ci_combined <- HPDinterval(mcmc_list, prob = 0.95)

# Convert the result to a data frame for easier plotting
ci_df <- as.data.frame(ci_combined)
ci_df$Parameter <- rownames(ci_df)
ci_df <- ci_df[, c("lower", "upper", "Parameter")]
colnames(ci_df) <- c("Lower", "Upper", "Parameter")

# Calculate the 95% credible intervals for resampled particles (ASMC algorithms)
ci_resampled1 <- apply(resampled_particles1, 2, function(x) quantile(x, c(0.025, 0.975)))
ci_resampled2 <- apply(resampled_particles2, 2, function(x) quantile(x, c(0.025, 0.975)))

ci_resampled_df1 <- data.frame(Parameter = colnames(resampled_particles1),
                               Lower = ci_resampled1[1, ],
                               Upper = ci_resampled1[2, ])

ci_resampled_df2 <- data.frame(Parameter = colnames(resampled_particles2),
                               Lower = ci_resampled2[1, ],
                               Upper = ci_resampled2[2, ])

# Add Algorithm column to each data frame
ci_df$Algorithm <- "AMCMC"
ci_resampled_df1$Algorithm <- "ASMC-AMTM"
ci_resampled_df2$Algorithm <- "AASMC"

# Combine all data frames into one
ci_all <- rbind(ci_df, ci_resampled_df1, ci_resampled_df2)

true_values <- c(β = 0.2976, δ = 1/3.2, γ = 1/8.5, ϕ = 100, S0 = 986, E0 = 8, I0 = 6)

library(ggplot2)

# Assuming ci_all has the correct structure
# Rename columns if necessary
colnames(ci_all) <- c("Lower", "Upper", "Parameter", "Algorithm")

# Update the parameter names if you have specific names for them
parameter_names <- c("β", "δ", "γ", "ϕ", "S0", "E0", "I0")
ci_all$Parameter <- factor(ci_all$Parameter, levels = parameter_names)

# Create a data frame for true values
true_values_df <- data.frame(Parameter = names(true_values), TrueValue = true_values)
true_values_df$Parameter <- factor(true_values_df$Parameter, levels = parameter_names)

# Create the plot
p <- ggplot(ci_all, aes(x = Algorithm, ymin = Lower, ymax = Upper, color = Algorithm)) +
  geom_pointrange(aes(y = (Lower + Upper) / 2), position = position_dodge(width = 0.5)) +
  facet_wrap(~Parameter, scales = "free_y") +
  geom_hline(data = true_values_df, aes(yintercept = TrueValue), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(title = "",
       y = "Parameter Value",
       x = NULL,
       color = "Method") +
  theme(axis.text.x = element_blank(),  # Remove x-axis text labels
        axis.ticks.x = element_blank(), # Remove x-axis ticks
        legend.position = "bottom",
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 14),  # Increase legend text font size
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  scale_color_manual(values = c("red", "green", "blue"))

# Save the plot as a high-resolution JPEG
ggsave("CI_SEIR.jpg", plot = p, width = 12, height = 8, dpi = 300)
