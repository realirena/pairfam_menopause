
rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores(logical= FALSE))
model_dir <- "U:/Documents/repos/pairfam_menopause/"
results_dir <- "U:/Documents/repos/pairfam_menopause/samples/"
seed = 1028

dt <- read.csv("meno_affect_06172026.csv")
traj_dt <- read.csv("affect_traj_06172026.csv")

N = nrow(traj_dt)
I = length(unique(dt$new_id))

other_cov = as.matrix(traj_dt[,c("nkids", "yeduc", "mar_stat2", "ethnic", "ever_smk")])
out_cov = as.matrix(dt[,c("baseline_kids","baseline_ed", "baseline_marstat2", "ethnic", "ever_smk")])
n_cov = ncol(other_cov)
K= 2

y = as.matrix(traj_dt[,c(1:10)])
P = ncol(y)


##  from the model 
compiled_model <- stan_model(paste0(model_dir, "lfm_joint_2k.stan"))

model_out <- sampling(compiled_model,
                      # include = TRUE,
                      sample_file=paste0(results_dir, '2lf_doublecov_recode_0705.csv'), #writes the samples to CSV file
                      iter =4000,
                      warmup=2000, #BURN IN
                      save_warmup=FALSE,
                      chains =4,
                      seed = seed,
                      control = list(max_treedepth = 60,
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
                        Tstart = dt$age1,
                        Tend = dt$age2,
                        n_cov = n_cov,
                        other_cov = other_cov,
                        out_cov = out_cov
                      ))

