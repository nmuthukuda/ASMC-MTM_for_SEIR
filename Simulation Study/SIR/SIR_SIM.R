set.seed(444)

# Solving I(t) DE trajectory using ODE solver
library(deSolve)
library(mvtnorm)
library(MASS)
library(coda)
library(MCMCpack)
library(truncnorm)
library(tmvtnorm)

# initial (state) values for SIR model
N <- 1000
x0 <- c(S = N-6, I = 6, R = 0)

# vector of time steps
times <- 0:59

# vector of parameters used in SIR model
params <- c(beta = 0.2976, gamma = 1/8.5)

SIR <- function(t, x, params) {
  with(as.list(c(params, x)), {
    dS <- -beta * S * I / N 
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I 
    list(c(dS, dI, dR))
  })
}

r1 <-  rk4(x0, times, SIR, params)

I_values <- r1[, "I"]

(I_max <- max(I_values))
mu_t <- I_values

h <- 100
p1 <- h/(h + mu_t)
data <- rnbinom(60, size = h, mu = mu_t)

data_table <- data.frame(Time = times, I_Values = round(I_values),Y = data)

# Convert the result to a data frame
r1_df <- as.data.frame(r1)

# Extract I values and calculate mu_t
I_values <- r1_df$I
mu_t <- I_values

# Generate observed data
h <- 100
p1 <- h / (h + mu_t)
data <- rnbinom(60, size = h, prob = p1)

# Create a data frame for the observed data
data_table <- data.frame(Time = times, I_Values = round(I_values), Y = data)

# Create a data frame for the SIR model with a variable to distinguish each component
r1_long <- r1_df %>%
  pivot_longer(cols = c(S, I, R), names_to = "Compartment", values_to = "Population")

# Set the factor levels for Compartment to ensure the desired order in the legend
r1_long$Compartment <- factor(r1_long$Compartment, levels = c("S", "I", "R"))

# Plot the SIR model trajectories and observations using ggplot2 with increased linewidth and legend
cplot1 <- ggplot() +
  geom_line(data = r1_long, aes(x = time, y = Population, color = Compartment), size = 1) +
  geom_point(data = data_table, aes(x = Time, y = Y, shape = "Observations"), color = "black", size = 1) +
  scale_color_manual(values = c(S = "blue", I = "red", R = "green"), 
                     labels = c(S = "S", I = "I", R = "R")) +
  scale_shape_manual(values = c(Observations = 16), labels = c(Observations = "Observations")) +
  labs(title = "", x = "Time", y = "Population", color = "", shape = "") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) +
  guides(shape = guide_legend(order = 2), color = guide_legend(order = 1))

# Save the plot
ggsave("SIR_SIM.jpeg", plot = cplot1, width = 12, height = 8, dpi = 300)

cplot1