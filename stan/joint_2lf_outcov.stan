//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  matrix get_lambda_matrix(vector lambda, int K, int P, int pos_loads){
    matrix[P,K] Lambda  = rep_matrix(0, P, K); //initiate matrix to zeros 
    Lambda[1, 1] = 1;
    Lambda[2:pos_loads,1] = lambda[2:pos_loads]; //positive loadings 
    Lambda[pos_loads+1,2] = 1;
    Lambda[pos_loads+2:P,2] = lambda[pos_loads+2:P]; // negative loadings -> this is set up to work on 2 latent factors
    return Lambda;
  }
  
  matrix get_rest_diag_matrix(matrix Sigma, int K){
    matrix[K,K] Sigma0 = Sigma; 
    for(k in 1:K){
     Sigma0[k,k] = 1; 
    }
    return(Sigma0);
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // total sample size 
  int P; // no. of observed outcomes 
  int K; // no. of latent factors 
  int<lower=1> I; //no. of individuals 
  int pos_loads; // keeps track of which emotions load to the positive affect 
  int<lower=1> ids[N]; //keeps track of which observations belong to which individuals 
  matrix[N, P] y;
  vector[N] time; // basis expansion on time (for the latent trajectories)
  // outcome model
  vector<lower=0>[I] Tstart;
  vector<lower=0>[I] Tend;
  vector[I] cens; // indicator if censored or not 
  int<lower=1> n_cov; // number of extra covariates for the survival regression
  matrix[I, n_cov] out_cov; 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector<lower=0.00001>[P] lambda; // factor loadings vectors that map the latent trajectories to the observed data
  real<lower=0> sigma_lambda;
  vector[K] B; 
  corr_matrix[2] Omega_k[K];  
  vector<lower=0>[2] tau_k[K]; 
  matrix[2,K] ran_eff[I]; // random effects for the latent trajectories 
  vector<lower=0>[P] S;
  // outcome parameters
  matrix[K, 2] b_rf; 
  vector[n_cov] phi_out; 
  real<lower=0> c;
  real b0; 
  
}

transformed parameters {
  matrix[P,K] Lambda;
  matrix[2,2] Sigma_k[K]; 
  vector[I] hazard;
  Lambda = get_lambda_matrix(lambda, K, P, pos_loads);

  for(k in 1:K){
    Sigma_k[k] =  diag_pre_multiply(tau_k[k],Omega_k[k]) * diag_pre_multiply(tau_k[k], Omega_k[k])';
  }

  for(i in 1:I){
   hazard[i] = exp(-(b0 + sum(ran_eff[i] .* b_rf) + dot_product(out_cov[i,],phi_out)) * c); // 
  } 
 // print(Lambda);
}

// The model to be estimated. 

model {
  matrix[N,P] y_mu;
  matrix[N, K] eta; // latent trajectories
  
  //weibull
  c ~ lognormal(1, 1);
  b0 ~ normal(0,1); 
  phi_out ~ normal(0,1);
  //latent factor
  lambda ~ normal(1, sigma_lambda);

  
  sigma_lambda ~ cauchy(0,2.5);
  S ~ cauchy(0, 2.5);
  
  //regression coefficients
  for(k in 1:K){
    to_vector(b_rf[k]) ~ normal(0,1);
    B[k] ~ normal(0, 5);
    Omega_k[k] ~ lkj_corr(1);
    tau_k[k] ~ cauchy(0,1);
  }
  
  ## have negative affect be drawn from normal w/ positive mean 
  //random effects
  for(i in 1:I){
    ran_eff[i][,1] ~ multi_normal(rep_vector(0,2), Sigma_k[1]); 
    ran_eff[i][,2] ~ multi_normal(rep_vector(0,2), Sigma_k[2]); 
  }
  
  
  for(n in 1:N){
    for(k in 1:K){
      eta[n,k] =B[k]*time[n] +  ran_eff[ids[n]][1,k] + ran_eff[ids[n]][2,k]*time[n]; //evolution of the latent trajectories  + dot_product(phi[k,], other_cov[ids[n],])
    }
  }
  y_mu = eta*Lambda';
  for(n in 1:N){
    for(p in 1:P){
        y[n,p] ~ normal(y_mu[n,p], S[p]);
    }
  }
 // outcome model
  for (i in 1:I) {
    if(cens[i]==1){
      target += weibull_lccdf(Tend[i] | c, hazard[i]); // if right censored
     //print(log_diff_exp(weibull_lcdf(Tend[i] | c, lambda), weibull_lcdf(Tstart[i] | c, lambda)));
    } else if(cens[i] == 0) {
       // interval censored 
       target += log_diff_exp(weibull_lcdf(Tend[i] | c, hazard[i]),weibull_lcdf(Tstart[i] | c, hazard[i]));
    } else if(cens[i]==2){ //left censoring
      target += weibull_lcdf(Tend[i] |c, hazard[i]);
    }
  }

}

## predicted times for each individual in the dataset (can plot observed vs predicted)
generated quantities{
  vector[I] pred;
  vector[I] surv;

  vector[1000] post_pred_pa;
  vector[1000] post_pred_na;
  vector[1000] post_pred_pbo;
  
  real hazard_pa;
  real hazard_na;
  real hazard_pbo;
  real hazard_ratio_pa;
  real hazard_ratio_na;

    // generate hazard ratio for positive affect slopes
  hazard_pa = exp(-(b0 + b_rf[2, 2])*c);
  hazard_na = exp(-(b0 + b_rf[2,1])*c); 
  hazard_pbo = exp(-b0*c);
  hazard_ratio_pa = exp((b_rf[1,2])*c);
  hazard_ratio_na = exp(b_rf[1,1]*c);

  // generate survival times (for plotting survival curves)
  for(i in 1:1000){
    post_pred_pa[i] = weibull_rng(c,  hazard_pa);
    post_pred_na[i] = weibull_rng(c, hazard_na);
    post_pred_pbo[i] = weibull_rng(c,  hazard_pbo);
  }
  
  for(i in 1:I){
    pred[i] = weibull_rng(c,hazard[i]);
    if(cens[i] == 1) {
      // Right censored
      surv[i] = weibull_lccdf(Tend[i] | c, hazard[i]);
    } else if(cens[i] == 0) {
      // Interval censored
      surv[i] = log_diff_exp(weibull_lcdf(Tend[i] | c, hazard[i]),
                                 weibull_lcdf(Tstart[i] | c, hazard[i]));
    } else if(cens[i] == 2) {
      // Left censored
      surv[i] = weibull_lcdf(Tend[i] | c, hazard[i]);
    }
  }
}


