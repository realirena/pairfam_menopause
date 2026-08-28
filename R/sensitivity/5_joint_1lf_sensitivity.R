
rm(list=ls())
library(haven)
library(rstan)
options(mc.cores = parallel::detectCores(logical= FALSE))
model_dir <- "U:/Documents/repos/pairfam_menopause/R/sensitivity/"
data_dir <- "U:/Documents/repos/menopause_models/R/data/sensitivity/"
results_dir <- "G:/irena/lfm/samples/sensitivity/long_covs/"
seed = 1028

dt <- read.csv(paste0(data_dir, "meno_srh_07272026.csv"))
traj_dt <- read.csv(paste0(data_dir, "srh_traj_07272026.csv"))

N = nrow(traj_dt)
I = length(unique(dt$new_id))

other_cov = as.matrix(traj_dt[,c("nkids", "yeduc", "mar_stat", "ethnic", "ever_smk")])
out_cov = as.matrix(dt[,c("baseline_kids","baseline_ed", "baseline_marstat", "ethnic", "ever_smk")])
n_cov = ncol(other_cov)
P = 1

y = traj_dt$srh

compiled_model <- stan_model(paste0(model_dir, "joint_1lf_long_covs.stan"))

model_out <- sampling(compiled_model,
                      # include = TRUE,
                      sample_file=paste0(results_dir, 'joint_1lf_long_covs.csv'), #writes the samples to CSV file
                      iter =2000,
                      warmup=1000, #BURN IN
                      chains = 4,
                      seed = seed,
                      control = list(max_treedepth = 15,
                                     adapt_delta=0.99),
                      data = list(
                        N = N,
                        P = P,
                        I = I,
                        time = traj_dt$age_std,
                        pos_loads = 5,
                        y = y,
                        ids = traj_dt$new_id,
                        cens = dt$cens_type,
                        Tstart =dt$age1_s30,
                        Tend = dt$age2_s30,
                        n_cov = n_cov,
                        other_cov = other_cov,
                        out_cov = out_cov
                      ))
