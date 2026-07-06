# ──────────────────────────────────────────────────────────────────────────────
# 1. SETUP & LIBRARIES
# ──────────────────────────────────────────────────────────────────────────────
rm(list = ls())
set.seed(112124) # Reproducibility

# tidyverse loads ggplot2, dplyr, tidyr, purrr, readr, etc.
library(tidyverse)
library(rstan)
library(haven)
library(lubridate)
library(ggnewscale)
library(gridExtra) # Was used but not loaded

options(mc.cores = parallel::detectCores(logical = FALSE))

# ──────────────────────────────────────────────────────────────────────────────
# 2. PATHS & DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────
data_dir   <- "U:/Documents/repos/menopause_models/R/"
results_dir <- "G:/irena/lfm/samples/"

# Helper for cross-platform safe paths
data_file <- function(name) file.path(data_dir, "data", name)

meno_affect_df <- read.csv(data_file("meno_affect_06172026.csv"))
meno_srh_df    <- read.csv(data_file("meno_srh_06172026.csv"))
affect_df      <- read.csv(data_file("affect_traj_06172026.csv"))
srh_df         <- read.csv(data_file("srh_traj_06172026.csv"))

# ──────────────────────────────────────────────────────────────────────────────
# 3. STAN MODEL LOADING
# ──────────────────────────────────────────────────────────────────────────────
# Kept only the second vector (first was overwritten)
sample_file_names <- c("joint_1lf_doublecov_0617_4", "joint_1lf_doublecov_0617_3",
                       "joint_1lf_doublecov_0617_2", "joint_1lf_doublecov_0617_1")

model_out <- read_stan_csv(file.path(results_dir, paste0(sample_file_names, ".csv")))

# Separate draws (for plotting) from summaries (for medians/stats)
post_draws   <- extract(model_out)
ran_eff_sum  <- summary(model_out, pars = "ran_eff", probs = c(0.025, 0.975))$summary

# ──────────────────────────────────────────────────────────────────────────────
# 4. WEIBULL HELPERS & PLOTTING
# ──────────────────────────────────────────────────────────────────────────────
weibull_hazard   <- function(t, shape, scale) (shape / scale) * (t / scale)^(shape - 1)
weibull_survival <- function(t, shape, scale) exp(-(t / scale)^shape)

make_weibull_df <- function(shape_draws, scale_draws, group_label, grid, type = "hazard") {
  fn <- if (type == "hazard") weibull_hazard else weibull_survival
  map_dfr(seq_along(shape_draws), ~ {
    data.frame(draw = .x, group = group_label, time = grid,
               value = fn(grid, shape_draws[.x], scale_draws[.x]))
  })
}

time_grid <- seq(30, 60)

shared_theme <- theme_bw() +
  theme(
    plot.title    = element_text(size = 25, face = "bold",  hjust = 0.5),
    plot.subtitle = element_text(size = 20, face = "bold",  hjust = 0.5),
    axis.text.x   = element_text(size = 15, angle = 45, vjust = 1, hjust = 1),
    axis.text.y   = element_text(size = 15),
    axis.title    = element_text(size = 20),
    legend.text   = element_text(size = 15),
    legend.title  = element_text(size = 18),
    strip.text    = element_text(size = 14)
  )

colors_1factor <- c("SRH Intercept" = "#3e821b", "Baseline" = "#F08030")
colors_2factor <- c("Positive Affect Slope" = "#3e821b", "Baseline" = "#F08030")

weibull_plot <- function(df, y_var, y_label, color_values, facet = FALSE) {
  p <- ggplot(df, aes(x = time, y = .data[[y_var]],
                      group = interaction(group, draw), color = group)) +
    geom_line(alpha = 0.05, linewidth = 0.3) +
    scale_color_manual(values = color_values) +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) +
    labs(x = "Age", y = y_label, color = "Group") + shared_theme
  
  if (facet) p <- p + facet_wrap(~group) + theme(legend.position = "none")
  p
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. 1-FACTOR MODEL PLOTS
# ──────────────────────────────────────────────────────────────────────────────
haz_df_1f <- bind_rows(
  make_weibull_df(post_draws$c, post_draws$hazard_int, "SRH Intercept", time_grid, "hazard"),
  make_weibull_df(post_draws$c, post_draws$hazard_pbo, "Baseline",      time_grid, "hazard")
)
surv_df_1f <- bind_rows(
  make_weibull_df(post_draws$c, post_draws$hazard_int, "SRH Intercept", time_grid, "survival"),
  make_weibull_df(post_draws$c, post_draws$hazard_pbo, "Baseline",      time_grid, "survival")
)

p_haz_1f  <- weibull_plot(haz_df_1f,  "value", "Hazard of Entering Perimenopause", colors_1factor, facet = TRUE)
p_surv_1f <- weibull_plot(surv_df_1f, "value", "Probability of Not Yet Entering Perimenopause", colors_1factor)

# ──────────────────────────────────────────────────────────────────────────────
# 6. 2-FACTOR MODEL PLOTS
# ──────────────────────────────────────────────────────────────────────────────
haz_df_2f <- bind_rows(
  make_weibull_df(post_draws$c, post_draws$hazard_pa,  "Positive Affect Slope", time_grid, "hazard"),
  make_weibull_df(post_draws$c, post_draws$hazard_pbo, "Baseline",              time_grid, "hazard")
)
surv_df_2f <- bind_rows(
  make_weibull_df(post_draws$c, post_draws$hazard_pa,  "Positive Affect Slope", time_grid, "survival"),
  make_weibull_df(post_draws$c, post_draws$hazard_pbo, "Baseline",              time_grid, "survival")
)

p_haz_2f  <- weibull_plot(haz_df_2f,  "value", "Hazard of Entering Perimenopause", colors_2factor, facet = TRUE)
p_surv_2f <- weibull_plot(surv_df_2f, "value", "Probability of Not Yet Entering Perimenopause", colors_2factor)

# ──────────────────────────────────────────────────────────────────────────────
# 7. MEDIAN & DIFFERENCE CALCULATIONS (DRY principle)
# ──────────────────────────────────────────────────────────────────────────────
compute_weibull_medians <- function(h1, h2, c, label1, label2) {
  med1 <- median(h1 * (-log(0.5))^(1/c))
  med2 <- median(h2 * (-log(0.5))^(1/c))
  diff_draws <- h2 * (-log(0.5))^(1/c) - h1 * (-log(0.5))^(1/c)
  
  cat(sprintf("Median age (%s): %.2f\n", label1, med1))
  cat(sprintf("Median age (%s): %.2f\n", label2, med2))
  cat(sprintf("Difference: %.2f\n", med2 - med1))
  cat(sprintf("95%% CrI: [%.2f, %.2f]\n\n", 
              quantile(diff_draws, 0.025), quantile(diff_draws, 0.975)))
  
  list(med1 = med1, med2 = med2, diff = med2 - med1, 
       cr_interval = quantile(diff_draws, c(0.025, 0.975)))
}

compute_weibull_medians(post_draws$hazard_pbo, post_draws$hazard_int, post_draws$c, 
                        "Baseline", "SRH Intercept")
compute_weibull_medians(post_draws$hazard_pbo, post_draws$hazard_pa,  post_draws$c, 
                        "Baseline", "PA Slope")

# ──────────────────────────────────────────────────────────────────────────────
# 8. INDIVIDUAL TRAJECTORY ANALYSIS
# ──────────────────────────────────────────────────────────────────────────────
# Extract & reshape random effects
ran_eff_df <- as.data.frame(ran_eff_sum)
ran_eff_df$affect   <- c("neg", "pos")
ran_eff_df$rf       <- rep(c("int", "slope"), each = 2)
ran_eff_df$affect_rf <- paste0(ran_eff_df$affect, "_", ran_eff_df$rf)
ran_eff_df$id       <- rep(meno_affect_df$new_id, each = 4)
ran_eff_df$censored_type <- rep(meno_affect_df$censored_type, each = 4)

plot_data <- ran_eff_df %>%
  select(id, mean, affect, rf, censored_type) %>%
  pivot_wider(names_from = rf, values_from = mean)

# Join with predicted ages
compare_data <- data.frame(summary(model_out, pars = "pred"))$summary
compare_data$id <- meno_affect_df$new_id

estimates <- left_join(
  compare_data[, c("id", "mean")],
  plot_data %>% filter(affect == "pos") %>% select(id, affect, int, slope),
  by = "id"
)

est_traj <- merge(estimates, affect_df[, c("new_id", "wave", "age_std", "age")],
                  by.x = "id", by.y = "new_id")

# Corrected trajectory calculation (removed duplicate assignment)
est_traj$pos_indiv_traj <- -0.06453948 + est_traj$int + (est_traj$slope * est_traj$age_std)

# Final plot (x = age, not wave)
c1 <- ggplot(data = est_traj, aes(x = age, y = pos_indiv_traj, group = as.factor(id))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.25) +
  geom_vline(xintercept = 41.7, linewidth = 1, linetype = "dashed", 
             alpha = 0.6, color = "gray40") +
  labs(x = "Age", y = "Standardized Positive Affect (Z-score)", color = "") +
  annotate("text", x = 40.78, y = 0.45, size = 14, 
           label = "Est. age at perimenopause: 46.78") + 
  annotate("text", x = 42.71, y = -0.75, size = 14, 
           label = "Est. age at perimenopause: 48.71") +
  shared_theme

# ──────────────────────────────────────────────────────────────────────────────
# 9. FINAL LAYOUT
# ──────────────────────────────────────────────────────────────────────────────
gridExtra::grid.arrange(p_surv_2f, c1, ncol = 1)