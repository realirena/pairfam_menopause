
// modify for 1 latent factor 
functions {
  matrix get_lambda_matrix(vector lambda, int K, int P, int pos_loads){
    matrix[P,K] Lambda  = rep_matrix(0, P, K); //initiate matrix to zeros
    Lambda[1, 1] = 1;
    Lambda[2:pos_loads,1] = lambda[2:pos_loads]; //positive loadings
   // Lambda[pos_loads+1,2] = 1;
   // Lambda[pos_loads+2:P,2] = lambda[pos_loads+2:P]; // negative loadings -> this is set up to work on 2 latent factors
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
  int<lower=1> I; //no. of individuals 
  int n_cov; // no. of other covariates
  int<lower=1> ids[N]; //keeps track of which observations belong to which individuals 
  vector[N] y;
  vector[N] time; // basis expansion on time (for the latent trajectories)
  matrix[N, n_cov] other_cov; 
//survival submodel
  vector<lower=0>[I] Tstart;
  vector<lower=0>[I] Tend;
  vector[I] cens; // indicator if censored or not 


}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> S;
  real<lower=0.00001> lambda; // factor loadings vectors that map the latent trajectories to the observed data
  real<lower=0> sigma_lambda;
  corr_matrix[2] Omega_k; 
  vector[2] ran_eff[I]; // random effects for the latent trajectories 
  real B; //population coefficient for the latent factor 
  vector<lower=0>[2] tau_k; //variance of the random effects 
  vector[n_cov] phi;
  // outcome parameters
  vector[2] b_rf; 
  real<lower=0> c;
  real b0; 
  
}

transformed parameters {
  matrix[2,2] Sigma_k;
  vector[I] hazard;
  
  Sigma_k  =  diag_pre_multiply(tau_k,Omega_k) * diag_pre_multiply(tau_k, Omega_k)';
  for(i in 1:I){
   hazard[i] = exp(-(b0 + sum(ran_eff[i] .* b_rf)) * c);
  } 
}



model {
  vector[N] eta; 
  vector[N] y_mu; 
  Omega_k ~ lkj_corr(1);
  tau_k ~ cauchy(0,1); 
  S ~ cauchy(0, 2.5);
  sigma_lambda ~ cauchy(0,1); 
  c ~ lognormal(1,1);
  lambda ~ normal(1, sigma_lambda);
  b0 ~ normal(0,1)
  B ~ normal(0, 1);
  phi ~ normal(0,1);
  for(i in 1:I){
    ran_eff[i] ~ multi_normal(rep_vector(0, 2), Sigma_k);
  }
  
  for(n in 1:N){
      y_mu[n] = B*time[n]  + dot_product(other_cov[n,], phi) + ran_eff[ids[n]][1] + ran_eff[ids[n]][2]*time[n]; // 
  }
  for(n in 1:N){
      y[n] ~ normal(y_mu[n], S);
  }
  
  // outcome model
  for (i in 1:I) {
    if(cens[i]==1){
      target += weibull_lccdf(Tend[i] | c, hazard[i]); // if right censored
     //print(log_diff_exp(weibull_lcdf(Tend[i] | c, lambda), weibull_lcdf(Tstart[i] | c, lambda)));
    } else if(cens[i] == 0) {
       // print(weibull_lccdf(Tend[i] | c, lambda));
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
 vector[1000] post_pred_int;
  vector[1000] post_pred_slope;
  vector[1000] post_pred_pbo;
  real hazard_int;
  real hazard_slope;
  real hazard_pbo;
  real hazard_ratio_int;
  real hazard_ratio_slope;
  
    // generate hazard ratio for positive affect slope
  hazard_int= exp(-(b0 + b_rf[1])*c);
  hazard_slope = exp(-(b0 + b_rf[2])*c);
  hazard_pbo = exp(-b0*c);
  hazard_ratio_int = exp((b_rf[1])*c);
  hazard_ratio_slope = exp((b_rf[2])*c);
  // generate survival times (for plotting survival curves)
  for(i in 1:1000){
    post_pred_int[i] = weibull_rng(c,  hazard_int);
    post_pred_slope[i] = weibull_rng(c,  hazard_slope);
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






