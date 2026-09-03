# =============================================================================
# Data preparation: PAIRFAM perimenopause joint models
#
# Produces, for each of the two models (SRH and depressiveness/affect):
#   - a person-level survival dataset  (one row per woman)
#   - a person-observation trajectory dataset (one row per woman-wave)
#
# Includes time-origin-shifted survival times (see ORIGINS below), so the Stan
# models can be run from a non-birth time origin without re-prepping the data.
#
# NOTE: several behavioural questions in the original script are flagged with
# "REVIEW:" below. Behaviour is preserved as-is; nothing was silently changed.
# =============================================================================

rm(list = ls())
library(survival)
library(tidyverse)
library(haven)

seed <- 112124
options(mc.cores = parallel::detectCores(logical = FALSE))

# ---- configuration ----------------------------------------------------------

DATA_DIR   <- "./data"
OUT_DIR    <- "./data/sensitivity"
REPO_DIR   <- "U:/Documents/repos/menopause_models/R/data/sensitivity"
DATE_TAG   <- format(Sys.Date(), "%m%d%Y")

# Time origins to include as shifted columns (age1_s30/age2_s30, age1_s35/...).
# 30 is the primary specification; 35 is the sensitivity check.
ORIGINS <- c(30, 35)

# Negative values in PAIRFAM encode different kinds of missingness.
neg_to_na <- function(x) replace(x, which(x < 0), NA)

# ---- 1. survival (person-level) data ----------------------------------------

df          <- readRDS(file.path(DATA_DIR, "pairfam_upto14_cohort3women.rds"))
meno_smoke  <- readRDS("./scripts/pairfam_meno_interval.rds")
meno_raw    <- readRDS(file.path(DATA_DIR, "pairfam_meno_interval_days.rds"))

# 22 women have no smoking data and are dropped here.
meno_df <- meno_raw |>
  left_join(meno_smoke[, c("id", "ever_smoke")], by = "id") |>
  filter(!is.na(ever_smoke)) |>
  distinct() |>
  mutate(
    age1 = as.numeric(age1),
    age2 = as.numeric(age2),
    age1 = if_else(censored_type == "left_censored", age2, age1),
    int_length = age2 - age1
  )

# ---- 2. longitudinal (person-observation) data ------------------------------
anchor_vars <- c("hlt1", "hlt7", "id", "wave", "age", "nkids", "yeduc", "marstat")

# Depressiveness items are the per2* columns; hlt1 is SRH.
health_df <- df |>
  filter(id %in% meno_df$id) |>
  dplyr::select(contains("per2"), all_of(anchor_vars))

per2_cols <- names(health_df)[str_detect(names(health_df), "^per2")]
stopifnot(length(per2_cols) == 10)   # 5 negative + 5 positive affect items

# Ethnicity is fixed at first observation.
meno_ethni <- df |>
  filter(id %in% meno_df$id) |>
  dplyr::select(id, ethni) |>
  distinct() |>
  mutate(across(everything(), neg_to_na))

# Depressiveness starts at wave 2, but SRH is available from wave 1, so this
# discards a wave of usable SRH data. Consider applying `wave > 1` only to the
# affect subset if that was not intended.

health_traj <- health_df |>
  left_join(meno_df, by = "id") |>
  left_join(meno_ethni, by = "id") |>
  filter(age <= round(as.numeric(age1), 0)) |>         
  mutate(across(everything(), neg_to_na))

names(health_traj)[match(per2_cols, names(health_traj))] <- c(
  "melancholy", "happy", "depressed", "sad", "desperate",
  "gloomy", "good", "secure", "calm", "enjoy"
)
health_traj <- health_traj |> rename(srh = hlt1)
health_traj <- health_traj |> mutate(age = as.numeric(age))
# Order items so the 5 negative-affect items precede the 5 positive-affect
# items. This ordering is what `pos_loads = 5` in the Stan models refers to.
health_traj <- health_traj |> relocate(happy, .after = gloomy)

affect_label_vars <- c("melancholy", "happy", "depressed", "sad", "desperate",
                       "gloomy", "good", "secure", "calm", "enjoy")

# Covariate recoding: reference categories are unmarried/divorced/widowed,
# ethnic German (incl. part-German), and never-smoker.
health_traj <- health_traj |>
  mutate(
    mar_stat = if_else(marstat == 2, 1, 0),
    ethnic   = if_else(ethni %in% c(1, 2, 3), 0, 1),
    ever_smk = if_else(ever_smoke == 2, 0, 1)
  )


# ---- 3. build the two analytic samples --------------------------------------

shared_covs <- c("age", "nkids", "yeduc", "ever_smoke", "marstat", "ethni")
 
# One function replaces the four near-identical blocks in the original script:
# complete-case filtering, truncation at onset, new_id assignment, and the
# matching person-level survival dataset.
#
# lag_years: buffer between the last retained health observation and age1.
#   0 (default) = the main specification: observations up to onset.
#   3           = sensitivity check: observations must end >= 3 years before
#                 onset, so trajectories cannot reflect early perimenopausal
#                 symptoms. Addresses the reverse-causation concern raised in
#                 the Discussion (health measures possibly being an early
#                 manifestation of perimenopause rather than preceding it).
# min_obs: minimum usable observations per woman (2 = enough to inform a slope).
build_sample <- function(traj, outcome_vars, min_wave = 1, min_obs = 2,
                         lag_years = 0) {
 
  traj_sub <- traj |>
    filter(wave >= min_wave) |>
    filter(if_all(all_of(c(outcome_vars, shared_covs)), ~ !is.na(.x))) |>
    filter(age < age1 - lag_years) |>
    group_by(id) |>
    filter(n() >= min_obs) |>
    ungroup() |>
    arrange(id, wave)
 
  id_map <- tibble(id = unique(traj_sub$id), new_id = seq_along(unique(traj_sub$id)))
 
  # Standardize age within this sample; keep the constants, which are needed to
  # back-transform for figures (e.g. the age at which intercepts are evaluated).
  age_mean <- mean(as.numeric(traj_sub$age))
  age_sd   <- sd(as.numeric(traj_sub$age))
 
  traj_sub <- traj_sub |>
    left_join(id_map, by = "id") |>
    mutate(
      age     = as.numeric(age),
      age_std = (age - age_mean) / age_sd
    ) |>
    arrange(new_id, wave, age)

    # baseline covariates = each woman's first (earliest-age) observation,
  # matching the manuscript's "computed on women's baseline measurements"
  baseline_cov <- traj_sub |>
    group_by(new_id) |>
    slice_min(age, n = 1, with_ties = FALSE) |>
    ungroup() |>
    dplyr::select(new_id,
                  baseline_kids    = nkids,
                  baseline_ed      = yeduc,
                  baseline_marstat = mar_stat,
                  ethnic,
                  ever_smk)
  
  surv_sub <- meno_df |>
    inner_join(id_map, by = "id") |>
    left_join(baseline_cov, by = "new_id") |>
    arrange(new_id)
 
  list(traj = traj_sub, surv = surv_sub,
       age_mean = age_mean, age_sd = age_sd, lag_years = lag_years)
}
 
# Depressiveness items are only collected from wave 2 onward; SRH from wave 1.
affect <- build_sample(health_traj, affect_label_vars, min_wave = 2)
srh    <- build_sample(health_traj, "srh",             min_wave = 1)

# Sensitivity: trajectories must end at least 3 years before onset. (check for reverse causality)
affect_lag3 <- build_sample(health_traj, affect_label_vars, min_wave = 2, lag_years = 3)
srh_lag3    <- build_sample(health_traj, "srh",             min_wave = 1, lag_years = 3)

# ---- 4. censoring codes and time-origin shifts ------------------------------

# Stan convention: 0 = interval, 1 = right, 2 = left.
code_censoring <- function(x) {
  out <- case_when(
    str_detect(x, "^left_censored")  ~ 2,
    str_detect(x, "^right_censored") ~ 1,
    str_detect(x, "^interval")       ~ 0,
    TRUE ~ NA_real_
  )
  stopifnot(!any(is.na(out)))   # fail loudly on an unrecognised label
  out
}

add_survival_cols <- function(surv, origins = ORIGINS) {
  out <- surv |> mutate(cens_type = code_censoring(censored_type))

  for (o in origins) {
    out[[paste0("age1_s", o)]] <- out$age1 - o
    out[[paste0("age2_s", o)]] <- out$age2 - o
  }
  out
}
affect <- build_sample(health_traj, affect_label_vars, min_wave = 2)
srh    <- build_sample(health_traj, "srh", min_wave = 1)

affect$surv <- add_survival_cols(affect$surv)
srh$surv    <- add_survival_cols(srh$surv)

affect_lag3$surv <- add_survival_cols(affect_lag3$surv)
srh_lag3$surv  <- add_survival_cols(srh_lag3$surv)
# > c(srh$age_mean, srh$age_sd)
# + c(affect$age_mean, affect$age_sd)
# [1] 41.153838  3.258025
# [1] 41.733482  2.961656


table(health_traj$wave)   # wave 1 should appear with ~966
nrow(srh$traj)            # up from 8,778
nrow(srh$surv)            # watch for a change from 932

# ---- 5. validation ----------------------------------------------------------

check_sample <- function(s, label) {
  cat("\n---", label, "---\n")
  cat("women:", nrow(s$surv), " observations:", nrow(s$traj), "\n")
  cat("censoring (0=int, 1=right, 2=left):\n")
  print(table(s$surv$cens_type))
 
  # new_id must run 1..n and match between the two files, or the Stan model
  # will pair women's survival times with the wrong random effects.
  stopifnot(
    identical(sort(unique(s$traj$new_id)), seq_len(nrow(s$surv))),
    !any(is.na(s$surv$age2)),
    all(s$surv$age1 <= s$surv$age2)
  )
 
  # Weibull requires strictly positive times on every shifted scale.
  for (o in ORIGINS) {
    a1 <- s$surv[[paste0("age1_s", o)]]
    a2 <- s$surv[[paste0("age2_s", o)]]
    cat("origin", o, "- min age1_s:", round(min(a1, na.rm = TRUE), 2),
        " min age2_s:", round(min(a2, na.rm = TRUE), 2), "\n")
    if (min(c(a1, a2), na.rm = TRUE) <= 0) {
      warning("Origin ", o, " produces non-positive times in ", label)
    }
  }
 
  cat("NAs in trajectory data:\n")
  print(colSums(is.na(s$traj))[colSums(is.na(s$traj)) > 0])
  invisible(NULL)
}
 
check_sample(affect, "affect / depressiveness")
check_sample(srh,    "SRH")
 
cat("\nIDs in SRH but not affect:", setdiff(srh$traj$id, affect$traj$id), "\n")
 
# ---- 6. write out -----------------------------------------------------------

# Writes both a copy retaining `id` (local, for traceability) and a copy with
# `id` dropped (repo, for analysis).
write_pair <- function(d, stem) {
  write_csv(d, file.path(OUT_DIR, paste0(stem, "_", DATE_TAG, ".csv")))
  write_csv(dplyr::select(d, -id),
            file.path(REPO_DIR, paste0(stem, "_", DATE_TAG, ".csv")))
}

write_pair(affect$surv, "meno_affect")
write_pair(srh$surv,    "meno_srh")
write_pair(affect$traj, "affect_traj")
write_pair(srh$traj,    "srh_traj")



# -----------------------------------------------------------
### pull baseline values for Susie
# ---- observed values at first observation ----
baseline_srh <- srh$traj |>
  group_by(new_id, id) |>
  slice_min(age, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(id, new_id, wave, age, srh)

baseline_affect <- affect$traj |>
  group_by(new_id, id) |>
  slice_min(age, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(id, new_id, wave, age,
         all_of(affect_label_vars))     # the ten items

# a summed/averaged score may be more useful than ten separate items
baseline_affect <- baseline_affect |>
  mutate(
    neg_affect_mean = rowMeans(across(c(melancholy, depressed, sad, desperate, gloomy))),
    pos_affect_mean = rowMeans(across(c(happy, good, secure, calm, enjoy)))
  )

write.csv(baseline_affect, paste0("./data/baseline/baseline_affect.csv"), row.names=FALSE)
write.csv(baseline_srh, paste0("./data/baseline/baseline_srh.csv"), row.names=FALSE)

# =============================================================================
# Optional: bounded left-censoring intervals
#
# Treats left-censored cases as interval-censored on (a_min, age2] rather than
# (0, age2], where a_min is the earliest plausible onset age. Produces extra
# columns rather than overwriting, so both specifications stay available.
# =============================================================================

add_bounded_cols <- function(surv, a_min_age = 40, origin = 30) {
  a_min <- a_min_age - origin
  age1_s <- surv[[paste0("age1_s", origin)]]

  surv[[paste0("Tstart_b", origin)]] <- if_else(surv$cens_type == 2, a_min, age1_s)
  surv[[paste0("cens_b")]]           <- if_else(surv$cens_type == 2, 0, surv$cens_type)
  surv
}


# affect$surv <- add_bounded_cols(affect$surv, a_min_age = 40, origin = 30)
# srh$surv    <- add_bounded_cols(srh$surv,    a_min_age = 40, origin = 30)

### turnbull 
srh$surv <- srh$surv |>
  mutate(
    cens_type = code_censoring(censored_type),
    tb_lower  = if_else(cens_type == 1, age2, if_else(cens_type == 0, age1, 0)),
    tb_upper  = if_else(cens_type == 1, NA_real_, age2)
  )

fit_tb <- survfit(Surv(tb_lower, tb_upper, type = "interval2") ~ 1, data = srh$surv)

plot_dat <- data.frame(
  time  = fit_tb$time,
  surv  = fit_tb$surv,
  lower = fit_tb$lower,
  upper = fit_tb$upper
)
plot_dat$cumhaz       <- -log(plot_dat$surv)
plot_dat$cumhaz_lower <- -log(plot_dat$upper)
plot_dat$cumhaz_upper <- -log(plot_dat$lower)



tail(plot_dat[, c("time", "surv", "lower", "upper")], 10)
min(plot_dat$surv)
plot(fit_tb, xlab = "Age", ylab = "Survival", xlim = c(40, 55))

plot_dat <- plot_dat[is.finite(plot_dat$cumhaz) & !is.na(plot_dat$lower), ]
range(plot_dat$time)
max(plot_dat$cumhaz)

t0_candidate <- min(plot_dat$time[plot_dat$surv < 1])
head(plot_dat[, c("time", "surv")], 10)

# corrected Turnbull, shifted scale
plot_dat_s <- transform(plot_dat, time = time - origin)
t0_s <- 44 - origin


post <- rstan::extract(model_out)
quantile(as.numeric(post$c), c(0.025, 0.5, 0.975))
ncol(post$hazard)

ORIGIN <- 30
age_grid   <- seq(44, 49.8, by = 0.1)
age_grid_s <- age_grid - ORIGIN
t0_s       <- 44 - ORIGIN

# corrected Turnbull, shifted and renormalized to age 44
tb_s <- transform(plot_dat, time = time - ORIGIN)
S_44 <- approx(tb_s$time, tb_s$surv, xout = t0_s)$y
tb_s$cumhaz       <- -log(tb_s$surv  / S_44)
tb_s$cumhaz_lower <- -log(tb_s$upper / S_44)
tb_s$cumhaz_upper <- -log(tb_s$lower / S_44)
tb_s <- tb_s[tb_s$time >= t0_s & is.finite(tb_s$cumhaz), ]

draws_w <- list(hazard = post$hazard, c = as.numeric(post$c))

res_corrected <- compare_model_to_turnbull(weibull_surv_fn, draws_w, age_grid_s, tb_s, t0 = t0_s)
res_corrected$summary

cf <- res_corrected$model_dat
cf$turnbull_cumhaz <- approx(tb_s$time, tb_s$cumhaz, xout = cf$time)$y
cf$surv_diff <- exp(-cf$cumhaz) - exp(-cf$turnbull_cumhaz)
cf$age <- cf$time + ORIGIN

summary(cf$surv_diff)
max(abs(cf$surv_diff), na.rm = TRUE)

plot(cf$age, cf$surv_diff, type = "b", xlab = "Age",
     ylab = "Model minus Turnbull (survival scale)")
abline(h = 0, lty = 2)


