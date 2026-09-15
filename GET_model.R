library(ggplot2)
library(reshape2) 

# Simulation settings
t <- 20160
dt <- 0.01



gen_suicide_sim <- function(t, dt, initial_params, other_params) {

  ## starting values
  S0 = initial_params[1]
  A0 = initial_params[2] 
  U0 = initial_params[3]
  T0 = initial_params[4] 
  O0 = initial_params[5] 
  E0 = initial_params[6] 
  I0 = initial_params[7]
  
  # params
  sigma = other_params["sigma"]
  f1 = other_params["f1"] 
  K2 = other_params["K2"]  
  b2 = other_params["b2"]  
  a2 = other_params["a2"]  
  d2 = other_params["d2"]  
  e2 = other_params["e2"]  
  g2 = other_params["g2"]  
  c3 = other_params["c3"]  
  b3 = other_params["b3"]  
  d4 = other_params["d4"]  
  c41 = other_params["c41"]  
  c42 = other_params["c42"]  
  e5 = other_params["e5"]  
  c51 = other_params["c51"]  
  c52 = other_params["c52"]  
  K6 = other_params["K6"]  
  f6 = other_params["f6"]  
  b6 = other_params["b6"]  
  c6 = other_params["c6"]  
  K7 = other_params["K7"]  
  g7 = other_params["g7"] 
  b7 = other_params["b7"]  
  c7 = other_params["c7"] 
  
  mu = (sigma^2) / 2
  
  # creating empty variable vectors to be filled up later, together with the starting values
  stressor <- rep(0, t); stressor[1] = S0
  av_state <- rep(0, t); av_state[1] = A0
  urge_escape <- rep(0, t); urge_escape[1] = U0
  sui_thoughts <- rep(0, t); sui_thoughts[1] = T0
  other_escape <- rep(0, t); other_escape[1] = O0
  ext_change <- rep(0, t); ext_change[1] = E0
  int_change <- rep(0, t); int_change[1] = I0
  
  
  for (i in 1:(t-1)){
    stressor[i + 1] = max(stressor[i] * exp((mu - (sigma^2)/2) * dt + sigma * rnorm (1,0, sqrt(dt)) - f1 * ext_change[i]), 0)
    av_state[i + 1] = max(av_state[i] + dt * (b2 * av_state[i] * (K2 - av_state[i]) + a2 * stressor[i] - d2 * sui_thoughts[i] - e2 * other_escape[i] - g2 * int_change[i]), 0)
    urge_escape[i + 1] = max(urge_escape[i] + dt * (-c3 * urge_escape[i] + b3 * av_state[i]), 0)
    sui_thoughts[i + 1] = max(sui_thoughts[i] + dt * (-d4 * sui_thoughts[i] + (1 / (1 + exp(-c41 * (urge_escape[i] - c42))))), 0)
    other_escape[i + 1] = max(other_escape[i] + dt * (-e5 * other_escape[i] + (1 / (1 + exp(-c51 * (urge_escape[i] - c52))))), 0)
    ext_change[i + 1] = max(ext_change[i] + dt * (f6 * ext_change[i] * (K6 - ext_change[i]) + b6 * av_state[i] - c6 * urge_escape[i]), 0)
    int_change[i + 1] = max(int_change[i] + dt * (g7 * int_change[i] * (K7 - int_change[i]) + b7 * av_state[i] - c7 * urge_escape[i]), 0)
  }
  
  
  
  output <- cbind(stressor, av_state, urge_escape, sui_thoughts, other_escape, ext_change, int_change)
  colnames(output) <- c("S", "A", "U", "T", "X", "E", "I")
  return(output)
}


# starting values
initial_params <- c(0.2, 0.3, 0.25, 0, 0.2, 0.05, 0.1)

# model parameters 
other_params <- c(sigma=0.1, f1=0.0001, K2=0.2, b2=4, a2=2, d2=1.5, e2=1, g2=0.5,
                  c3=3, b3=1.5, d4=1, c41=100, c42=0.25, e5=3, c51=50, c52=0.2,
                  K6=0.1, f6=0.5, b6=0.41, c6=0.82, K7=0.05, g7=0.5, b7=0.65, c7=1.3)

# running the model
set.seed(504) 
sim <- gen_suicide_sim(t, dt, initial_params, other_params)

# reshape data for plot
sim_df <- as.data.frame(sim)
sim_df$time <- seq(0, (t-1)*dt, by = dt)
sim_long <- melt(sim_df, id.vars = "time", variable.name = "Variable", value.name = "Value")

# Plot
ggplot(sim_long, aes(x = time, y = Value, color = Variable)) +
  geom_line(linewidth = 0.8) +
  labs(title = "1 Iteration of the model",
       x = "Time",
       y = "Value") +
  theme_bw() +
  theme(text = element_text(size = 14)) 




