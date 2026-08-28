
rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores(logical= FALSE))
model_dir <- "U:/Documents/repos/pairfam_menopause/R/sensitivity/"
data_dir <- "U:/Documents/repos/menopause_models/R/data/sensitivity/"
results_dir <- "G:/irena/lfm/samples/sensitivity/long_covs/"
seed = 1028

dt <- read.csv(paste0(data_dir, "meno_affect_07272026.csv"))
traj_dt <- read.csv(paste0(data_dir, "affect_traj_07272026.csv"))

N = nrow(traj_dt)
I = length(unique(dt$new_id))

other_cov = as.matrix(traj_dt[,c("nkids", "yeduc", "mar_stat", "ethnic", "ever_smk")])
out_cov = as.matrix(dt[,c("baseline_kids","baseline_ed", "baseline_marstat", "ethnic", "ever_smk")])
n_cov = ncol(other_cov)
K= 2

y = as.matrix(traj_dt[,c(1:10)])
P = ncol(y)


##  from the model 
compiled_model <- stan_model(paste0(model_dir, "joint_2lf_long_covs.stan"))

model_out <- sampling(compiled_model,
                      # include = TRUE,
                      sample_file=paste0(results_dir, '2lf_long_covs.csv'), #writes the samples to CSV file
                      iter =2000,
                      warmup=1000, #BURN IN
                      save_warmup=FALSE,
                      chains =4,
                      seed = seed,
                      control = list(max_treedepth = 12,
                                     adapt_delta=0.95),
                      data = list(
                        I = I,
                        N = N,
                        K = K,
                        P = P,
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

