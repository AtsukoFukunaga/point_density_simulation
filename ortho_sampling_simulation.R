library(tidyverse)

load("point_density/point_density_data.rda")  # specify appropriate path to the point_density_data.rda file

# simulation code

# create vectors to simulate the number of photoquads and the number of annotation points
n_photoquads_hi <- c(10, 20, 40, 60, 80, 100, 120)  # Hawaii Island
n_photoquads_kapou <- c(10, 20, 30, 40)  # Kapou
annotation_points <- c(10, 20, 50, 100, 200, 400, 600, 800, 1000)

# create functions for bootstrapping
# option 1: get a point estimate for each photoquad and use the estimate for bootstrapping
# option 2: combine all the annotation points from photoquads for bootstrapping

# function for bootstrapping - option 1
boot_stat_quad <- function(boot_dat, seed) {
  
  boot_quad <- boot_dat %>% 
    group_by(name) %>% 
    summarise(live_coral_quad = mean(live_coral))
  mean_live_coral_quad <- mean(boot_quad$live_coral_quad)
  se_live_coral_quad <- sd(boot_quad$live_coral_quad) / sqrt(nrow(boot_quad))
  n <- nrow(boot_quad)
  bt_df <- data.frame(rep = numeric(), stat = numeric(), t = numeric())
  set.seed(seed)
  for (l in 1:10000) {
    vec <- 1:n
    boot_vec <- sample(vec, n, replace = TRUE)
    boot_live_coral_quad <- boot_quad$live_coral_quad[boot_vec]
    boot_live_coral_quad_mean <- mean(boot_live_coral_quad)
    boot_live_coral_quad_se <- sd(boot_live_coral_quad) / sqrt(length(boot_live_coral_quad))
    boot_live_coral_quad_t <- (boot_live_coral_quad_mean - mean_live_coral_quad) / boot_live_coral_quad_se
    temp_bt_df <- data.frame(rep = l, stat = boot_live_coral_quad_mean, t = boot_live_coral_quad_t)
    bt_df <- bind_rows(bt_df, temp_bt_df)
  }
  boot_stat_quad_mean <- mean(bt_df$stat)
  boot_t_quad_qt <- quantile(bt_df$t, probs = c(0.025, 0.975))
  boot_quad_ci_l <- mean_live_coral_quad - boot_t_quad_qt[[2]] * se_live_coral_quad
  boot_quad_ci_u <- mean_live_coral_quad - boot_t_quad_qt[[1]] * se_live_coral_quad
  boot_stat_quad_summary <- c(boot_stat_quad_mean, boot_quad_ci_l, boot_quad_ci_u)
  
  return(boot_stat_quad_summary)
}

# function for bootstrapping - option 2
boot_stat_plot <- function(boot_dat, seed) {
  
  mean_live_coral <- mean(boot_dat$live_coral)
  se_live_coral <- sd(boot_dat$live_coral) / sqrt(nrow(boot_dat))
  n <- nrow(boot_dat)
  bt_df <- data.frame(rep = numeric(), stat = numeric(), t = numeric())
  set.seed(seed)
  for (l in 1:10000) {
    vec <- 1:n
    boot_vec <- sample(vec, n, replace = TRUE)
    boot_live_coral <- boot_dat$live_coral[boot_vec]
    boot_live_coral_mean <- mean(boot_live_coral)
    boot_live_coral_se <- sd(boot_live_coral) / sqrt(length(boot_live_coral))
    boot_live_coral_t <- (boot_live_coral_mean - mean_live_coral) / boot_live_coral_se
    temp_bt_df <- data.frame(rep = l, stat = boot_live_coral_mean, t = boot_live_coral_t)
    bt_df <- bind_rows(bt_df, temp_bt_df)
  }
  boot_stat_plot_mean <- mean(bt_df$stat)
  boot_t_qt <- quantile(bt_df$t, probs = c(0.025, 0.975))
  boot_ci_l <- mean_live_coral - boot_t_qt[[2]] * se_live_coral
  boot_ci_u <- mean_live_coral - boot_t_qt[[1]] * se_live_coral
  boot_stat_plot_summary <- c(boot_stat_plot_mean, boot_ci_l, boot_ci_u)
  
  return(boot_stat_plot_summary)
}

# subset data by plot

## West Hawaii plot names
## WH_2022_S11_169 (4.57%), WH_2022_S07_260 (9.01%), WH_2022_S04_020 (14.0%), 
## WH_2022_S02_026 (15.2%), WH_2022_S04_107 (20.4%), WH_2022_S04_155 (21.5%), 
## Kapou plot names
## LIS_CW1_T1_2021 (18.3%), LIS_10_T1_2021 (21.1%), LIS_COURT_T1_2021 (22.6%), 
## LIS_R10_T1_2021 (30.8%), LIS_R7_T1_2021 (38.4%)

plot_name <- "LIS_R7_T1_2021"  # enter plot name
n_photoquads <- n_photoquads_kapou  # enter n_photoquads_hi for Hawaii Island or n_photoquads_kapou for Kapou
simu_df <- simulation_data %>% 
  filter(plot == plot_name) %>% 
  select(plot, name, coral)
photoquad_names <- unique(simu_df$name)

# generate a dataset for simulation for the plot and perform bootstrapping
# repeat 1000 times with each one (round) having a unique random seed
# change the range of random seed numbers to ensure each round has a unique random seed
# running 1 round at a time - change "size" and adjust the range to run multiple rounds 
rand_seeds <- sample(1:10, size = 1, replace = FALSE)

# set up a dataframe to store results
boot_df <- data.frame(plot = character(0), 
                      n_photoquad = numeric(0), 
                      n_annotation_point = numeric(0), 
                      seed = numeric(0),
                      boot_plot_mean = numeric(0), 
                      boot_plot_025 = numeric(0), 
                      boot_plot_975 = numeric(0), 
                      boot_quad_mean = numeric(0), 
                      boot_quad_025 = numeric(0), 
                      boot_quad_975 = numeric(0))
counter <- 0

for (rand_seed in rand_seeds) {
  counter <- counter + 1
  print(counter)
  set.seed(rand_seed)  # set seed and make a note for reproducibility
  
  # generate a dataset
  simulated_plot <- list()
  
  for (i in 1:length(n_photoquads)) {
    n <- n_photoquads[i]
    simu_photoquad <- sample(photoquad_names, n, replace = FALSE)
    
    simulated_photoquads <- list()
    
    for (quad in simu_photoquad) {
      temp <- simu_df %>% filter(name == quad)
      
      simulated_points <- list()
      
      for (j in 1:length(annotation_points)) {
        pt <- annotation_points[j]
        vec_points <- 1:nrow(temp)
        simu_vec_points <- sample(vec_points, pt, replace = FALSE)
        temp2 <- temp[simu_vec_points, ]
        simulated_points[[j]] <- temp2
      }
      
      simulated_photoquads[[quad]] <- simulated_points
    }
    
    simulated_plot[[i]] <- simulated_photoquads
  }
  
  for (n_qd in 1:length(n_photoquads)) {
    print(n_photoquads[n_qd])
    sub_dat <- simulated_plot[[n_qd]]
    
    boot_dat <- data.frame(plot = character(0), name = character(0), coral = character(0))
    for (n_pt in 1:length(annotation_points)) {
      
      for (k in 1:length(sub_dat)) {
        boot_dat <- bind_rows(boot_dat, sub_dat[[k]][[n_pt]])
      }
      
      boot_dat <- boot_dat %>% 
        mutate(live_coral = ifelse(coral == "Other", 0, 1))
      
      boot_plot_sum <- boot_stat_plot(boot_dat, rand_seed)
      boot_quad_sum <- boot_stat_quad(boot_dat, rand_seed)
      
      temp_boot_df <- data.frame(plot = plot_name, 
                                 n_photoquad = n_photoquads[n_qd],  
                                 n_annotation_point = annotation_points[n_pt], 
                                 seed = rand_seed, 
                                 boot_plot_mean = boot_plot_sum[1], 
                                 boot_plot_025 = boot_plot_sum[2], 
                                 boot_plot_975 = boot_plot_sum[3], 
                                 boot_quad_mean = boot_quad_sum[1], 
                                 boot_quad_025 = boot_quad_sum[2], 
                                 boot_quad_975 = boot_quad_sum[3])
      
      boot_df <- bind_rows(boot_df, temp_boot_df)
    }
  }
}

# resulting boot_df dataframe contains:
# plot: plot name
# n_photoquad: the number of simulated photoquads
# n_annotation_point: the number of simulated annotation points per photoquad
# seed: unique seed
# boot_plot_mean: bootstrap mean coral cover under Option 2
# boot_plot_025: 0.025th quantile of bootstrapped samples under Option 2
# boot_plot_975: 0.975th quantile of bootstrapped samples under Option 2
# boot_quad_mean: bootstrap mean coral cover under Option 1
# boot_quad_025: 0.025th quantile of bootstrapped samples under Option 1
# boot_quad_975: 0.975th quantile of bootstrapped samples under Option 1

# save(boot_df, file = "point_density/pd_001.rda")
