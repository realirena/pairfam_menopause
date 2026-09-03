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
library(patchwork)
library(grid)
library(gridExtra) #

options(mc.cores = parallel::detectCores(logical = FALSE))

# ──────────────────────────────────────────────────────────────────────────────
# 2. PATHS & DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────
data_dir   <- "U:/Documents/repos/menopause_models/R"
results_dir <- "G:/irena/lfm/samples/sensitivity/lag3/"
ORIGIN = 30
# Helper for cross-platform safe paths
data_file <- function(name) file.path(data_dir, "data/sensitivity/", name)

meno_affect_df <- read.csv(data_file("meno_affect_07272026.csv"))
meno_srh_df    <- read.csv(data_file("meno_srh_07272026.csv"))
affect_df      <- read.csv(data_file("affect_traj_07272026.csv"))
srh_df         <- read.csv(data_file("srh_traj_07272026.csv"))

# ──────────────────────────────────────────────────────────────────────────────
# 3. STAN MODEL LOADING
# ──────────────────────────────────────────────────────────────────────────────
# Kept only the second vector (first was overwritten)
srh_stems    <- paste0("joint_1lf_0820_lag3_", 1:4)
affect_stems <- paste0("2lf_doublecov_origin35_", 1:4)

srh_model_out <- read_stan_csv(file.path(results_dir, paste0(srh_stems, ".csv")))
affect_model_out <- read_stan_csv(file.path(results_dir, paste0(affect_stems, ".csv")))

# ──────────────────────────────────────────────────────────────────────────────
# 3.1 estimate hazard ratios
# ──────────────────────────────────────────────────────────────────────────────
b   <- as.matrix(model_out, pars = "b_rf")     # draws x 4
tau <- as.matrix(model_out, pars = "tau_k")    # draws x 4
c_d <- as.numeric(as.matrix(model_out, pars = "c"))
colnames(b); colnames(tau)   # CHECK the ordering before pairing them

## one factor
# pair each coefficient with the tau for the same effect type and factor
std <- cbind(
  srh_int   = b[, "b_rf[1]"] * tau[, "tau_k[1]"],
  srh_slope   = b[, "b_rf[2]"] * tau[, "tau_k[2]"]
)

# two factor
# pair each coefficient with the tau for the same effect type and factor
std <- cbind(
  na_int   = b[, "b_rf[1,1]"] * tau[, "tau_k[1,1]"],
  pa_int   = b[, "b_rf[1,2]"] * tau[, "tau_k[2,1]"],
  na_slope = b[, "b_rf[2,1]"] * tau[, "tau_k[1,2]"],
  pa_slope = b[, "b_rf[2,2]"] * tau[, "tau_k[2,2]"]
)

# standardized coefficients (x100, as in the tables)
apply(std * 100, 2, quantile, c(0.025, 0.5, 0.975))

# time ratios -- note the c
apply(exp(-std * c_d), 2, quantile, c(0.025, 0.5, 0.975))

# what does origin 35 imply for the median age shift?
c_35   <- as.matrix(srh_model_out,, pars = "c")[, 1]
b0_35  <- as.matrix(srh_model_out,, pars = "b0")[, 1]
b_35   <- as.matrix(srh_model_out,, pars = "b_rf[1]")[, 1]
tau_35 <- as.matrix(srh_model_out,, pars = "tau_k[1]")[, 1]

med_ref <- exp(-b0_35 * c_35) * log(2)^(1/c_35) + 35
med_1sd <- exp(-(b0_35 + b_35 * tau_35) * c_35) * log(2)^(1/c_35) + 35

quantile(med_ref, c(0.025, 0.5, 0.975))          # expect ~46.1
quantile(med_1sd - med_ref, c(0.025, 0.5, 0.975)) # expect ~0.42

summary(srh_model_out, pars=c("c"))$summary
summary(affect_model_out, pars=c("c"))$summary


#coefficients trade off?
cor(as.matrix(srh_adjusted_out, pars = "b_rf"))[1,2]
cor(as.matrix(srh_longcov_out, pars = "b_rf"))[1,2]

# are the SRH random effects correlated with the survival covariates?
re_srh <- rstan::extract(srh_adjusted_out)$ran_eff
re_med <- apply(re_srh, c(2,3), median)
cor(re_med[,1], dt$baseline_ed); cor(re_med[,2], dt$baseline_ed)
cor(re_med[,1], dt$ever_smk);    cor(re_med[,2], dt$ever_smk)
cor(re_med[,2], dt$baseline_kids)
cor(re_med[,2], dt$baseline_marstat)
cor(re_med[,2], dt$ethnic)

obs_srh <- srh_df |>
  group_by(new_id) |>
  summarise(mean_srh = mean(srh), .groups = "drop")

int_df <- tibble(new_id = seq_along(re_med[,1]), int = re_med[,1]) |>
  left_join(obs_srh, by = "new_id")

cor(int_df$int, int_df$mean_srh)     # expect strongly positive
plot(int_df$int, int_df$mean_srh)

ed_df <- dt |> select(new_id, baseline_ed) |> left_join(obs_srh, by = "new_id")
cor(ed_df$baseline_ed, ed_df$mean_srh)

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

time_grid <- seq(1, 30)     # shifted scale, i.e. ages 31-55

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

colors_1factor <- c("SRH Intercept" = "#8653b2", "Baseline" = "#F08030")
colors_2factor <- c("Positive Affect Slope" = "#3e821b", "Baseline" = "#F08030")

weibull_plot <- function(df, y_var, y_label, color_values, facet = FALSE, origin_shift = 30) {
  p <- ggplot(df, aes(x = time + origin_shift, y = .data[[y_var]],
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

# hazard for a woman one SD above average on the SRH intercept
draws_mat <- as.matrix(srh_model_out, pars = c("b0", "b_rf", "tau_k", "c"))

b0_d    <- draws_mat[, "b0"]
b_int   <- draws_mat[, "b_rf[1]"]
tau_int <- draws_mat[, "tau_k[1]"]
c_d     <- draws_mat[, "c"]

hazard_pbo_calc <- exp(-b0_d * c_d)
hazard_int_calc <- exp(-(b0_d + b_int * tau_int) * c_d)

### sanity check
#diff_draws <- (hazard_int_calc - hazard_pbo_calc) * log(2)^(1/c_d)
#quantile(diff_draws, c(0.025, 0.5, 0.975))

#quantile(hazard_pbo_calc * log(2)^(1/c_d) + 30, c(0.025, 0.5, 0.975))   # expect ~46.12
#quantile(hazard_int_calc * log(2)^(1/c_d) + 30, c(0.025, 0.5, 0.975))   # expect ~46.55


haz_df_1f <- bind_rows(
  make_weibull_df(c_d, hazard_pbo_calc, "Baseline",      time_grid, "hazard"),
  make_weibull_df(c_d, hazard_int_calc, "SRH Intercept", time_grid, "hazard")
)

surv_df_1f <- bind_rows(
  make_weibull_df(c_d, hazard_pbo_calc, "Baseline",      time_grid, "survival"),
  make_weibull_df(c_d, hazard_int_calc, "SRH Intercept", time_grid, "survival")
)

p_haz_1f  <- weibull_plot(haz_df_1f,  "value", "Hazard of Entering Perimenopause", colors_1factor, facet = TRUE)
p_surv_1f <- weibull_plot(surv_df_1f, "value", "Probability of Not Yet Entering Perimenopause", colors_1factor)

# ──────────────────────────────────────────────────────────────────────────────
# 6. 2-FACTOR MODEL PLOTS
# ──────────────────────────────────────────────────────────────────────────────

dm <- as.matrix(affect_model_out, pars = c("b0", "b_rf", "tau_k", "c"))

b0_a   <- dm[, "b0"]
b_pa   <- dm[, "b_rf[2,2]"]
tau_pa <- dm[, "tau_k[2,2]"]
c_a    <- dm[, "c"]

hazard_pbo_a <- exp(-b0_a * c_a)
hazard_pa_a  <- exp(-(b0_a + b_pa * tau_pa) * c_a)


haz_df_2f <- bind_rows(
  make_weibull_df(c_a,  hazard_pa_a,  "Positive Affect Slope", time_grid, "hazard"),
  make_weibull_df(c_a,hazard_pbo_a, "Baseline",              time_grid, "hazard")
)
surv_df_2f <- bind_rows(
  make_weibull_df(c_a,  hazard_pa_a,"Positive Affect Slope", time_grid, "survival"),
  make_weibull_df(c_a, hazard_pbo_a, "Baseline",  time_grid, "survival")
)
## save as 800 x 350
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

# does the affect shift still come out at 0.61 with consistent extraction?
compute_weibull_medians(hazard_pbo_a, hazard_pa_a, c_a, "Baseline", "PA Slope")

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

# ──────────────────────────────────────────────────────────────────────────────
# 9a. alternative plot: plot medians
# ──────────────────────────────────────────────────────────────────────────────

make_weibull_summary <- function(shape_draws, scale_draws, group_label, grid,
                                 type = "survival", origin_shift = 30) {
  fn <- if (type == "hazard") weibull_hazard else weibull_survival

  # rows = draws, cols = grid points
  vals <- sapply(grid, function(t) fn(t, shape_draws, scale_draws))

  data.frame(
    group = group_label,
    time  = grid,
    age   = grid + origin_shift,
    med   = apply(vals, 2, median),
    lower = apply(vals, 2, quantile, 0.025),
    upper = apply(vals, 2, quantile, 0.975)
  )
}

time_grid <- seq(1, 30, by = 0.25)   # ages 31-60, finer grid for smooth curves

srh_surv_summary <- bind_rows(
  make_weibull_summary(c_d, hazard_pbo_calc, "Population average",     time_grid, "survival"),
  make_weibull_summary(c_d, hazard_int_calc, "One SD above average",   time_grid, "survival")
)

surv_summary <- bind_rows(
  make_weibull_summary(c_a, hazard_pbo_a, "Population average",       time_grid, "survival"),
  make_weibull_summary(c_a, hazard_pa_a,  "One SD above average",     time_grid, "survival")
)

shared_colors <- c("Population average" = "#F08030", "One SD above average" = "#3e821b")


g1 <- ggplot(srh_surv_summary, aes(x = age, y = med, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = shared_colors) +
  scale_fill_manual(values = shared_colors) +
  coord_cartesian(xlim = c(40, 56)) +
  labs(x = NULL, y = NULL, title = "Self-rated health",
       color = "Group", fill = "Group") +
  shared_theme

g2<- ggplot(surv_summary, aes(x = age, y = med, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = shared_colors) +
  scale_fill_manual(values = shared_colors) +
  coord_cartesian(xlim = c(40, 56)) +
  labs(x = "Age", y = "Probability of Not Yet Entering Perimenopause",title = "Depressiveness",
       color = "Group", fill = "Group") +
  shared_theme

g1 <- g1 + labs(y = "Survival probability")
g2 <- g2 + labs(y = "Survival probability")


(g1 / g2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.margin = margin(5, 5, 5, 20))


# ────────────────────────────────────────────────────────────────────────────── 
##### 10: reproduce figure 3 (comparison of two women)
# ──────────────────────────────────────────────────────────────────────────────
set.seed(12345)
 
# ---- 1. draws, kept on one extraction path so they stay aligned -------------
 
# rstan::extract() permutes draws by default; as.matrix() does not. Mixing the
# two silently misaligns iterations, so everything below uses as.matrix().
dm  <- as.matrix(affect_model_out, pars = c("b0", "b_rf", "tau_k", "c"))
hz  <- as.matrix(affect_model_out, pars = "hazard")
re  <- as.matrix(affect_model_out, pars = "ran_eff")
 
c_a       <- dm[, "c"]
tau_pa_sl <- dm[, "tau_k[2,2]"]   # factor 2 (positive affect), element 2 (slope)
 
stopifnot(identical(meno_affect_df$new_id, sort(meno_affect_df$new_id)))   # dt row order must match Stan's
 
# ---- 2. each woman's posterior median PA random effects ---------------------
 
# ran_eff is indexed [individual, effect type, factor]:
#   [i,1,2] = positive affect intercept, [i,2,2] = positive affect slope
I <- nrow(meno_affect_df)
 head(colnames(re), 8)
I
nrow(meno_affect_df)
pa_int   <- sapply(seq_len(I), function(i) median(re[, sprintf("ran_eff[%d,1,2]", i)]))
pa_slope <- sapply(seq_len(I), function(i) median(re[, sprintf("ran_eff[%d,2,2]", i)]))

# ---- 3. model-implied median onset age --------------------------------------

onset_age <- function(i, probs = c(0.025, 0.5, 0.975)) {
  quantile(hz[, sprintf("hazard[%d]", i)] * log(2)^(1 / c_a) + ORIGIN, probs)
}
 
# ---- 4. pick two comparable women -------------------------------------------
 
tau_med <- median(tau_pa_sl)
 
cand <- tibble(
  new_id   = seq_len(I),
  int      = pa_int,
  slope    = pa_slope,
  slope_sd = pa_slope / tau_med,          # slope in SD units
  n_obs    = as.integer(table(affect_df$new_id)[as.character(seq_len(I))])
)
 
# Aim for similar intercepts (near the median) but slopes about +1 and -1 SD.
# Requiring a reasonable number of observations avoids women whose random
# effects are mostly prior.
int_window <- quantile(cand$int, c(0.4, 0.6))
 
pick <- cand |>
  left_join(meno_affect_df |> select(new_id, baseline_kids, baseline_ed,
                         baseline_marstat, ethnic, ever_smk), by = "new_id") |>
  filter(between(int, int_window[1], int_window[2]), n_obs >= 8)

# among candidates, find pairs matching on the binary covariates
hi_pool <- pick |> filter(slope_sd > 0.7, slope_sd < 1.3)
lo_pool <- pick |> filter(slope_sd < -0.7, slope_sd > -1.3)
# then choose one from each with identical marstat / smoking / ethnicity
# 111, 73 match on 
woman_hi <- pick |> filter(new_id == 73)
woman_lo  <- pick |> filter(new_id == 111)

chosen <- bind_rows(woman_hi, woman_lo)
print(chosen)
 
# how far apart are they, and what onset does the model imply?
onset <- sapply(chosen$new_id, onset_age)
colnames(onset) <- chosen$new_id
print(round(onset, 2))
b0_a  <- dm[, "b0"]
b_pa  <- dm[, "b_rf[2,2]"]

# holding covariates at zero, varying only the PA slope random effect
med_at <- function(slope_val) {
  quantile(exp(-(b0_a + b_pa * slope_val) * c_a) * log(2)^(1/c_a) + 30,
           c(0.025, 0.5, 0.975))
}

# ---- 5. fitted trajectories over each woman's observed ages -----------------
 
est_traj <- affect_df |>
  filter(new_id %in% chosen$new_id) |>
  select(new_id, wave, age, age_std) |>
  left_join(chosen |> select(new_id, int, slope, slope_sd), by = "new_id") |>
  mutate(pa_traj = int + slope * age_std)
 
onset_q <- sapply(chosen$new_id, onset_age)
colnames(onset_q) <- chosen$new_id

labels <- tibble(
  new_id    = chosen$new_id,
  onset_med = onset_q["50%", ],
  lower     = onset_q["2.5%", ],
  upper     = onset_q["97.5%", ],
  label     = sprintf("Est. age at onset: %.1f (95%% CrI %.1f, %.1f)",
                      onset_med, lower, upper)
) |>
  left_join(
    est_traj |> group_by(new_id) |> summarise(x = max(age), y = last(pa_traj), .groups = "drop"),
    by = "new_id"
  )
# ---- 6. figure --------------------------------------------------------------
labels$x <- min(est_traj$age)
labels$y <- c(max(est_traj$pa_traj), min(est_traj$pa_traj))

ggplot(est_traj, aes(x = age, y = pa_traj,
                     group = factor(new_id), color = factor(new_id))) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_vline(xintercept = mean(affect_df$age), linetype = "dashed", color = "grey50") +
geom_text(data = labels, aes(x = x, y = y, label = label, color = factor(new_id)),
          hjust = -0.65, size = 6, show.legend = FALSE) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  labs(x = "Age", y = "Positive affect (standardized)", color = NULL) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        axis.text = element_text(size = 15))

### alternate without color: 
ggplot(est_traj, aes(x = age, y = pa_traj, group = factor(new_id))) +
  geom_line(aes(linetype = factor(new_id)), linewidth = 1.1, color = "black") +
  geom_point(size = 2, color = "black") +
  geom_vline(xintercept = mean(affect_df$age), linetype = "dashed", color = "grey50") +
  geom_text(data = labels, aes(x = x, y = y, label = label),
            hjust = -0.65, size = 6, color = "black") +
  scale_linetype_manual(values = c("solid", "dotted")) +
  labs(x = "Age", y = "Positive affect (standardized)") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        axis.text = element_text(size = 15))
### save as 1100 x 900

# =============================================================================
# Checks worth running before using the figure:
#
#   # 1. are the two women's intercepts actually similar?
#   chosen$int
#
#   # 2. how extreme are their slopes, in SD units? The original figure used
#   #    women at roughly +5 SD and 0 SD, which is not a representative contrast.
#   chosen$slope_sd
#
#   # 3. does the labelled onset differ from the posterior predictive mean that
#   #    the earlier version used?
#   sapply(chosen$new_id, function(i) mean(hz[, sprintf("hazard[%d]", i)] *
#            gamma(1 + 1/c_a)) + ORIGIN)
# =============================================================================
 