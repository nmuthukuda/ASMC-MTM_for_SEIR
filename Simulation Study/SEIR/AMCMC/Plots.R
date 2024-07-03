setwd("C://SFU//Research//Comparison//AdaptiveASMC//Results_Sim//SEIR//AMCMC")
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the CSV files
chain1 <- read.csv("AMCMC_Finale_chain1.csv")
chain2 <- read.csv("AMCMC_Finale_chain2.csv")
chain3 <- read.csv("AMCMC_Finale_chain3.csv")
chain4 <- read.csv("AMCMC_Finale_chain4.csv")
chain5 <- read.csv("AMCMC_Finale_chain5.csv")

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

gelman.diag(mcmc_list)

summary(mcmc_list)

effectiveSize(mcmc_list)


chain1 <- chain11%>% mutate(Iteration = 1:34000)
chain2 <- chain22%>% mutate(Iteration = 1:34000)
chain3 <- chain33%>% mutate(Iteration = 1:34000)
chain4 <- chain44%>% mutate(Iteration = 1:34000)
chain5 <- chain55%>% mutate(Iteration = 1:34000)

# Add a chain identifier column
chain1$Chain <- "Chain 1"
chain2$Chain <- "Chain 2"
chain3$Chain <- "Chain 3"
chain4$Chain <- "Chain 4"
chain5$Chain <- "Chain 5"

# Combine all chains into one data frame
combined_chains <- bind_rows(chain1, chain2, chain3, chain4, chain5)

# Rename the parameter columns
parameter_names <- c("beta", "delta", "gamma", "phi", "S0", "E0", "I0")
names(combined_chains)[1:7] <- parameter_names

# Reshape the data into a long format
long_chains <- combined_chains %>%
  pivot_longer(cols = all_of(parameter_names), names_to = "Parameter", values_to = "Value")

# Create a named vector for parsed labels
parsed_labels <- c(
  beta = expression(beta),
  gamma = expression(gamma),
  delta = expression(delta),
  phi = expression(phi),
  S0 = expression(S[0]),
  E0 = expression(E[0]),
  I0 = expression(I[0])
)

p <- ggplot(long_chains, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(linewidth = 0.8) +  # Increase the line width
  facet_wrap(~ Parameter, scales = "free_y", labeller = labeller(Parameter = label_parsed)) +
  labs(x = "Iteration", y = "Value") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),# Increase y-axis label font size
    axis.title.x = element_text(size = 16),  # Increase y-axis label font size
    legend.text = element_text(size = 14),  # Increase legend text font size
    legend.title = element_text(size = 16),  # Increase legend title font size
    panel.spacing = unit(1, "lines"),  # Adds space between facets
    strip.background = element_blank(),  # Remove background from facet labels
    strip.placement = "outside",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),  # Customize major grid lines
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),  # Customize minor grid lines
    axis.ticks = element_line(linewidth = 0.5),  # Customize axis ticks
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Add plot margins
  ) +
  guides(color = guide_legend(nrow = 1))


# Display the plot
print(p)

# Save the plot as a high-resolution JPEG
ggsave("trace_plots_SEIR.jpg", plot = p, width = 12, height = 8, dpi = 300)
