#### Load packages ####
library(sensobol)
library(deSolve)
library(doParallel)
library(foreach)
library(tidyverse)
library(progress)


### for Pawn 
library(SAFER)
library(progressr) 
library(caTools)
library(calibrater)
library(ggplot2)
library(gridExtra)
library(matrixStats)


# ___________________________________________________________________________________
######## Sobol ########-------------------------------------------------------------- 
# ___________________________________________________________________________________


###### Setup ### --------------------------------------------------------- 

N <- 2^11 # defining sample size


# using sobol estimator
R <- 10 ^ 3
type <- "norm"
conf <- 0.95


########### Sampling matrices ### ------------------------------------------------ 

params <- c(
  "sigma", "f1", "K2", "b2", "a2", "d2", "e2", "g2", 
  "c3", "b3", "d4", "c41", "c42", "e5", "c51", "c52", 
  "K6", "f6", "b6", "c6", "K7", "g7", "b7", "c7", "dummy"
)

mat <- sobol_matrices(
  matrices = c("A", "B", "AB"),
  N = N, 
  params = params, # parameters
  order = "first", # which order
  type = "QRN" # Quasi random numbers (sobol)
)

# Fix all parameters except d4, c42, c3, b3, a2, d2 to nominal values
mat[, "sigma"] <- 0.1
mat[, "f1"]    <- 0.0001
mat[, "K2"]    <- 0.2
mat[, "b2"]    <- 4
mat[, "a2"]    <- qunif(mat[, "a2"], 1, 10)
mat[, "d2"]    <- qunif(mat[, "d2"], 0.1, 0.5)
mat[, "e2"]    <- 1
mat[, "g2"]    <- 0.5
mat[, "c3"]    <- qunif(mat[, "c3"], 1, 5)
mat[, "b3"]    <- qunif(mat[, "b3"], 0.5, 2)
mat[, "d4"]    <- qunif(mat[, "d4"], 0.1, 5)
mat[, "c41"]   <- 100
mat[, "c42"]   <- qunif(mat[, "c42"], 0.1, 0.5)
mat[, "e5"]    <- 3
mat[, "c51"]   <- 50
mat[, "c52"]   <- 0.2
mat[, "K6"]    <- 0.1
mat[, "f6"]    <- 0.5
mat[, "b6"]    <- 0.41
mat[, "c6"]    <- 0.82
mat[, "K7"]    <- 0.05
mat[, "g7"]    <- 0.5
mat[, "b7"]    <- 0.65
mat[, "c7"]    <- 1.3
mat[, "dummy"] <- qunif(mat[, "dummy"], 0, 1)



##### Run  model ##### ----------------------------------------

# starting values to run the model
y0 = c(0.4, 0.2, 0.3, 0.1, 0.1, 0.1, 0.1) 


# Wrapper function to run model row-wise + progress bar
simulate_s1_batch_max_version <- function(param_matrix, t = t, dt = dt, y0 = initial_params) {
  
  # param_matrix: each row = a set of parameters
  n_runs <- nrow(param_matrix)
  results_scalar <- numeric(n_runs) 
  
  # progress bar
  pb <- progress_bar$new(
    format = "  Running simulations [:bar] :percent ETA: :eta",
    total = n_runs, 
    width = 50
  )
  
  # going through all the rows of param_matrix and calculate model
  for (i in 1:n_runs) {
    params <- as.numeric(param_matrix[i, ])
    names(params) <- colnames(param_matrix)
    
    # Run the simulation for this row of parameters
    set.seed(504) 
    sim_out <- gen_suicide_sim(t, dt, initial_params = y0, other_params = params)
    results_scalar[i] <- max(sim_out[, "T"]) 
    
    pb$tick()
  }
  
  return(results_scalar)
}

# run model over samples to get output
y <- simulate_s1_batch_max_version(mat, t = 2000, dt = 0.01, y0 = y0)
#saveRDS(y, "samples_sobol.rds")



###########  Output Distribution ### ------------------------------------------------ 

plot_uncertainty(Y = y, N = N)


########### Sobol Indices ### ------------------------------------------------ 

# parameters of interest
keep_params <- c("d4", "c42", "c3", "b3", "a2", "d2", "dummy")

# Calculate indices
ind <- sobol_indices(matrices = c("A", "B", "AB"), Y = y, N = N, params = params, boot = TRUE, R = R,
                     parallel = "no", type = type, conf = conf)

res <- ind$results %>% filter(parameters %in% keep_params)

light_colors <- c("Si" = "#7BAFD4",  # soft blue
                  "Ti" = "#A084C9")  # soft purple

res %>% 
  ggplot(aes(x = parameters, y = original, fill = sensitivity)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, color = "black") +   # bars
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci),
                width = 0.3,
                position = position_dodge(width = 0.7),
                color = "gray15") +           # errorbars
  scale_fill_manual(values = light_colors,
                    labels = c("Si" = "First-order",         # map original factor levels to new labels
                               "Ti" = "Total-order"),
                    #labels = c("First-order" = "Si", "Total-order" = "Ti"),
                    name = "Sobol index") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0,0)) +
  labs(x = "Parameter", y = "Sobol index") +
  ylim(0,1)+
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.position.inside = c(0.15, 0.78),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  )




########### Convergence ### ------------------------------------------------ 

sub.sample.sobol <- 2^(seq(8,11, 1)) # Sample Sizes for convergence

### sobol convergence 
cc <- sobol_convergence(matrices   = matrices <- c("A", "B", "AB"),
                        Y          = y,      
                        N          = N,
                        sub.sample = sub.sample.sobol,      
                        first      = "saltelli",
                        total      = "jansen",
                        order      = "first",
                        params     = params,
                        plot.order = T,
                        boot       = TRUE,
                        R          = 10)



##### Plot convergence ### ------------------------------------------------ 
par(mfrow= c(1,1))
### changing data for plot 
cc_plot_data <- cc$convergenc %>% 
  filter(parameters %in% keep_params) %>%
  mutate(sample_size = log2((Cost / (25 + 2))))
unique(cc_plot_data$Cost)

# Define parameters and sample sizes
keep_params <- c("a2", "d2", "c3", "b3", "d4", "c42", "dummy") # Adjust if needed
sample_sizes <- unique(cc_plot_data$Cost)# Adjust if Cost differs
xvals <-log2((sample_sizes / (25 + 2)))

# Function to prepare data for plotting
prepare_sobol_data <- function(data, sensitivity_type, params, sample_sizes) {
  # Filter data
  plot_data <- data[sensitivity == sensitivity_type & parameters %in% params]
  
  # Create matrices using tapply
  mean_matrix <- 
    tapply(plot_data$original, list(plot_data$Cost, plot_data$parameters), mean)
  
  lb_matrix <- 
    tapply(plot_data$low.ci, list(plot_data$Cost, plot_data$parameters), mean)
  
  
  ub_matrix <-
    tapply(plot_data$high.ci, list(plot_data$Cost, plot_data$parameters), mean)
  
  return(list(mean=mean_matrix, lb=lb_matrix, ub=ub_matrix))
}

# somwthing is going wrong here with the names
# Prepare data
ti_data <- prepare_sobol_data(cc$convergence, "Ti", keep_params, sample_sizes)
si_data <- prepare_sobol_data(cc$convergence, "Si", keep_params, sample_sizes)



# Plot settings from PAWN plot
colors <- 1:7 # Black, red, green, blue, cyan, magenta, yellow
pchs <- rep(15, 7) # Square, circle, triangle, etc.

# Plot for Total Sensitivity Indices (Ti)
par(mar = c(5, 4, 4, 6) + 0.1) # Match PAWN plot margins
plot(xvals, ti_data$mean[,1],
     type = "n",
     xlab = "Number of Samples",
     ylab = "Total Sensitivity Index",
     xlim = range(xvals),
     ylim = range(ti_data$lb, ti_data$ub, na.rm=TRUE),
     xaxt = "n")
axis(1, at = xvals, labels = parse(text = paste0("2^", xvals)))

# Ribbons
for (i in 1:ncol(ti_data$mean)) {
  polygon(c(xvals, rev(xvals)),
          c(ti_data$ub[,i], rev(ti_data$lb[,i])),
          col = adjustcolor(colors[i], alpha.f = 0.25), border = NA)
}

# Mean lines and points
matplot(xvals, ti_data$mean,
        type = "b", pch = pchs[1:ncol(ti_data$mean)], col = colors,
        lty = 1, lwd = 1.2, cex = 0.9, add = TRUE)

# Legend
par(xpd = NA)
legend("right",
       legend = colnames(ti_data$mean),
       col = colors, pch = pchs[1:ncol(ti_data$mean)], lty = 1,
       inset = c(-0.4, 0), bty = "n", cex = 0.9) # modify the -0.4 if legend is inside the plot







# ___________________________________________________________________________________
#### PAWN ########-------------------------------------------------------------- 
# ___________________________________________________________________________________


n <- 10 # number of conditioning intervals



##### Sampling matrix: ##### ----------------------------------------

# parameter range
DistrPar <- list(
  c(1, 10),           # a2      (random)
  c(0.1, 0.5),        # d2      (random)
  c(1, 5),            # c3      (random)
  c(0.5, 2),          # b3      (random)
  c(0.1, 5),          # d4      (random)
  c(0.1, 0.5),        # c42     (random)
  c(0, 1)             # dummy   (keep uniform for screening)
)

params <- c(
  "sigma", "f1", "K2", "b2", "a2", "d2", "e2", "g2", 
  "c3", "b3", "d4", "c41", "c42", "e5", "c51", "c52", 
  "K6", "f6", "b6", "c6", "K7", "g7", "b7", "c7", "dummy"
)

# the parameters that will be included
keep_params <- c("d4", "c42", "c3", "b3", "a2", "d2", "dummy")

# actual sampling matrix 
X <- AAT_sampling(samp_strat ="lhs" ,   # sampling strategy
                  M = length(keep_params),   # number of inputs 
                  distr_fun = "unif",      # distribution function
                  distr_par = DistrPar,    # parameter ranges
                  N = 2^14)               # number of samples
colnames(X) <- keep_params


##### Run  model ##### ----------------------------------------

# starting values for the model to run
y0 = c(0.4, 0.2, 0.3, 0.1, 0.1, 0.1, 0.1) # orig

# Wrapper function to run model rowwise + progress bar
simulate_s1_batch_max_version_PAWN_fixed <- function(param_matrix, t = t, dt = dt, y0 = initial_params) {
  
  # param_matrix: each row = a set of parameters
  n_runs <- nrow(param_matrix)
  results_scalar <- numeric(n_runs) 
  
  # progress bar
  pb <- progress_bar$new(
    format = "  Running simulations [:bar] :percent ETA: :eta",
    total = n_runs, 
    width = 50
  )
  
  # going through all the rows of param_matrix and calculate model
  for (i in 1:n_runs) {
    params <- as.numeric(param_matrix[i, ])
    names(params) <- colnames(param_matrix)
    
    # adding the fixed parameters into the model
    fixed <- c(sigma = 0.12,
               f1    = 0.0001,
               K2    = 0.2,
               b2    = 4,
               e2    = 1,
               g2    = 0.5,
               c41   = 100,
               e5    = 3,
               c51   = 50,
               c52   = 0.2,
               K6    = 0.1,
               f6    = 0.5,
               b6    = 0.41,
               c6    = 0.82,
               K7    = 0.05,
               g7    = 0.5,
               b7    = 0.65,
               c7    = 1.3)
    
    pars <- fixed
    pars["a2"]   <- param_matrix[i, "a2"]
    pars["d2"]   <- param_matrix[i, "d2"]
    pars["c3"]   <- param_matrix[i, "c3"]
    pars["b3"]   <- param_matrix[i, "b3"]
    pars["d4"]   <- param_matrix[i, "d4"]
    pars["c42"]  <- param_matrix[i, "c42"]
    pars["dummy"]<- param_matrix[i, "dummy"]
    
    
    # Run the simulation for this row of parameters
    set.seed(504) 
    sim_out <- gen_suicide_sim(t, dt, initial_params = y0, other_params = pars)
    results_scalar[i] <- max(sim_out[, "T"]) 
    
    pb$tick()
  }
  
  return(results_scalar)
}

# run model over samples to get output
y_pawn <- simulate_s1_batch_max_version_PAWN_fixed(X, t = 2000, dt = 0.01, y0 = y0)
#saveRDS(y_pawn, "samples_pawn.rds")


##### PAWN indices ##### ----------------------------------------

pawn_ind_new <- 
  with_progress(pawn_indices(X = X, 
                             Y = y_pawn, 
                             n = 10, 
                             Nboot = 20, 
                             dummy = T))




# computing the confidence intervlas
KS_stat <- pawn_ind_new$KS_mean  # extracting all the means
KS_max_m <- colMeans(KS_stat) #calculating the mean of all the means across the 10 intervals and the bootstap. 
alfa <- 0.05  # CI level
KS_max_lb <- colQuantiles(KS_stat,probs=alfa/2) # Lower bound of CI
KS_max_ub <- colQuantiles(KS_stat,probs=1-alfa/2) # Upper bound of CI


# put all the data together for the plot
dat_pawn_plot <- data.frame(keep_params,KS_max_m,KS_max_lb, KS_max_ub  )


##### Plot indices: ##### ----------------------------------------
ggplot(dat_pawn_plot, 
       aes(x = keep_params, 
           y = KS_max_m)) +
  geom_col(position = position_dodge(width = 1),
           width = 0.8, color = "black",
           fill = "#7BAFD4") +
  geom_errorbar(aes(ymin=KS_max_lb,
                    ymax=KS_max_ub, color = "black"), width = 0.6,
                position = position_dodge(width = 0.7),
                color = "gray1")+ 
  labs(x = "Parameter", y = "Sensitivity index")+
  theme_bw(base_size = 14)+  
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.78),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  )


##### Plot CDFs    ##### ----------------------------------------

M <- length(keep_params) # number of parameters, necessary to run the pawn_plot_CDF
pawn_cdf <- pawn_plot_CDF(X, y_pawn, n=10, labelinput=keep_params)
CDFs <- pawn_CDF(X, y_pawn, n = 10, dummy = T, verbose = TRUE)


##### Convergence ------------------------------------------------

N <- 2^16  # the original sample size
NN <- 2^(seq(8,11, 1)) # the sample sizes for convergence

# computing pawn indices for every sample size
pawn_conv <- pawn_convergence(X, y_pawn, n, NN, Nboot = 5)
KS_mean_c <- pawn_conv$KS_mean # extracting the means

# Compute mean and CI of the indices across the bootstrap samples:
alfa <- 0.05 # Significance level for CI for the convergence plot

KS_stat <- KS_mean_c
KS_mean_c_m <- t(sapply(KS_stat,colMeans)) # mean
KS_mean_c_lb <-  t(sapply(KS_stat,colQuantiles,probs=alfa/2)) # Lower bound
KS_mean_c_ub <- t(sapply(KS_stat,colQuantiles,probs=1-alfa/2)) # Upper bound

# changing the column names, necessary for the legend
colnames(KS_mean_c_m) <- keep_params



##### Plot Convergence  ------------------------------------------------------------------
xvals <- log2(NN) # the x values

# empty plot setup
par(mar = c(5, 4, 4, 6) + 0.1)
plot(xvals, KS_mean_c_m[,1],
     type = "n",
     xlab = "Number of Samples",
     ylab = "mean KS",
     xlim = range(xvals),
     ylim = range(KS_mean_c_lb, KS_mean_c_ub),
     xaxt = "n")  # suppress x-axis

# x-axis with 2^k labels
axis(1, at = xvals, labels = parse(text = paste0("2^", 8:11)))

# ribbons
for(i in 1:ncol(KS_mean_c_m)){
  polygon(c(xvals, rev(xvals)),
          c(KS_mean_c_ub[,i], rev(KS_mean_c_lb[,i])),
          col = adjustcolor(i, alpha.f = 0.25), border = NA)
}

# mean lines
matlines(xvals, KS_mean_c_m,
         type = "b", pch = 16, col = 1:7, lty = 1, lwd = 1.2, cex = 0.9)

# legend
par(xpd = NA)
legend("right",
       legend = colnames(KS_mean_c_m),
       col = 1:7, pch = 16, lty = 1,
       inset = c(-0.4, 0), bty = "n", cex = 0.9) # if legend is inside the plot, then adjust the -0.4


