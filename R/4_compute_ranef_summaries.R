# =============================================================================
# Per-woman random effect summaries
#
# Extracts posterior mean and 95% CrI for each woman's random intercept and
# slope, on both the raw latent scale and standardized by the corresponding
# random-effect SD (tau).
#
# NOTE ON EXTRACTION: everything comes from as.matrix(), which returns draws in
# chain order. rstan::extract() permutes them by default, and mixing the two
# silently misaligns iterations. Do not combine outputs from the two paths.
# =============================================================================

library(tidyverse)
library(rstan)

# ──────────────────────────────────────────────────────────────────────────────
# 2. PATHS & DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────
data_dir   <- "U:/Documents/repos/menopause_models/R"
results_dir <- "U:/Documents/repos/menopause_models/R/data/"
ORIGIN = 30
# Helper for cross-platform safe paths
data_file <- function(name) file.path(data_dir, "data/sensitivity/", name)
# ──────────────────────────────────────────────────────────────────────────────
# 3. STAN MODEL LOADING
# ──────────────────────────────────────────────────────────────────────────────
# Kept only the second vector (first was overwritten)
srh_stems <- paste0("joint_1lf_0720_origin_shift_left_", 1:4)
affect_stems    <- paste0("2lf_doublecov_0814_shift_left_", 1:4)

srh_model_out <- read_stan_csv(file.path(results_dir, paste0(srh_stems, ".csv")))
affect_model_out <- read_stan_csv(file.path(results_dir, paste0(affect_stems, ".csv")))


# ---- two-factor (depressiveness) model --------------------------------------
#
# ran_eff is indexed [individual, effect type, factor]:
#   [i,1,1] NA intercept   [i,1,2] PA intercept
#   [i,2,1] NA slope       [i,2,2] PA slope
#
# tau_k is indexed [factor, effect type] -- the transposed order:
#   tau_k[1,1] NA intercept SD   tau_k[1,2] NA slope SD
#   tau_k[2,1] PA intercept SD   tau_k[2,2] PA slope SD

summarise_ranef_2f <- function(fit, n_women) {

  re  <- as.matrix(fit, pars = "ran_eff")
  tau <- as.matrix(fit, pars = "tau_k")

  # (effect label, ran_eff index, matching tau column)
  spec <- tribble(
    ~effect,        ~re_idx,  ~tau_col,
    "na_intercept", "1,1",    "tau_k[1,1]",
    "pa_intercept", "1,2",    "tau_k[2,1]",
    "na_slope",     "2,1",    "tau_k[1,2]",
    "pa_slope",     "2,2",    "tau_k[2,2]"
  )

  map_dfr(seq_len(nrow(spec)), function(r) {
    idx <- spec$re_idx[r]
    tau_draws <- tau[, spec$tau_col[r]]

    map_dfr(seq_len(n_women), function(i) {
      draws     <- re[, sprintf("ran_eff[%d,%s]", i, idx)]
      std_draws <- draws / tau_draws     # per draw, so tau uncertainty propagates

      tibble(
        new_id   = i,
        effect   = spec$effect[r],
        mean     = mean(draws),
        lower    = unname(quantile(draws, 0.025)),
        upper    = unname(quantile(draws, 0.975)),
        mean_sd  = mean(std_draws),
        lower_sd = unname(quantile(std_draws, 0.025)),
        upper_sd = unname(quantile(std_draws, 0.975))
      )
    })
  })
}

# ---- one-factor (SRH) model -------------------------------------------------
#
# ran_eff is [individual, effect type]; tau_k is a plain length-2 vector.

summarise_ranef_1f <- function(fit, n_women) {

  re  <- as.matrix(fit, pars = "ran_eff")
  tau <- as.matrix(fit, pars = "tau_k")

  spec <- tribble(
    ~effect,         ~re_idx, ~tau_col,
    "srh_intercept", "1",     "tau_k[1]",
    "srh_slope",     "2",     "tau_k[2]"
  )

  map_dfr(seq_len(nrow(spec)), function(r) {
    tau_draws <- tau[, spec$tau_col[r]]

    map_dfr(seq_len(n_women), function(i) {
      draws     <- re[, sprintf("ran_eff[%d,%s]", i, spec$re_idx[r])]
      std_draws <- draws / tau_draws

      tibble(
        new_id   = i,
        effect   = spec$effect[r],
        mean     = mean(draws),
        lower    = unname(quantile(draws, 0.025)),
        upper    = unname(quantile(draws, 0.975)),
        mean_sd  = mean(std_draws),
        lower_sd = unname(quantile(std_draws, 0.025)),
        upper_sd = unname(quantile(std_draws, 0.975))
      )
    })
  })
}

# ---- run --------------------------------------------------------------------

ranef_affect <- summarise_ranef_2f(affect_model_out, n_women = nrow(meno_affect_df))
ranef_srh    <- summarise_ranef_1f(srh_model_out,    n_women = nrow(meno_srh_df))

# sanity checks: the standardized effects should have SD near 1 across women,
# and the raw SDs should be close to the estimated tau values
ranef_affect |>
  group_by(effect) |>
  summarise(sd_raw = sd(mean), sd_std = sd(mean_sd), .groups = "drop")

quantile(as.matrix(affect_model_out, pars = "tau_k")[, "tau_k[2,2]"],
         c(0.025, 0.5, 0.975))   # compare against sd_raw for pa_slope

# NB: sd(mean) across women will be SMALLER than tau, because posterior means
# are shrunk toward zero. That is expected, not an error.

# ---- wide format, one row per woman -----------------------------------------

ranef_affect_wide <- ranef_affect |>
  select(new_id, effect, mean, lower, upper, mean_sd) |>
  pivot_wider(names_from = effect,
              values_from = c(mean, lower, upper, mean_sd))

write_csv(ranef_affect, paste0(data_dir,  "/data/ranef_affect_summaries.csv"))
write_csv(ranef_srh,   paste0(data_dir, "/data/ranef_srh_summaries.csv"))

### join the data to the original Pairfam ids:
ranef_affect <- read.csv(paste0(data_dir, "/data/ranef_affect_summaries.csv"))
ranef_srh <- read.csv(paste0(data_dir, "/data/ranef_srh_summaries.csv"))

## affect IDs and SRH ids respectively (from the pairfam dataset)
meno_affect_df <-  read.csv("C:/cloud/Pairfam/data/sensitivity/meno_affect_07272026.csv")
meno_srh_df    <- read.csv("C:/cloud/Pairfam/data/sensitivity/meno_srh_07272026.csv")

affect_ids <- unique(meno_affect_df[,c("id", "new_id")])
srh_ids <- unique(meno_srh_df[,c("id", "new_id")])

ranef_affect <- merge(affect_ids, ranef_affect, by=c("new_id"))
srh_affect <- merge(srh_ids, ranef_srh, by=c("new_id"))

write_csv(ranef_affect, paste0(data_dir,  "/data/ranef_affect_summaries_original_id.csv"))
write_csv(srh_affect,   paste0(data_dir, "/data/ranef_srh_summaries_original_id.csv"))
