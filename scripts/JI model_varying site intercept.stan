data {
    int<lower=0> n_obs;
    int<lower=0> n_day;
    int<lower=0> n_site;
    int siteID[n_obs];
    vector[n_obs] X_pre; 
    vector[n_obs] X_t; 
    matrix[n_obs, n_day] W;
}
parameters {
    real b_0_mean;
    real<lower=0> b_0_sigma;
    vector[n_site] b_0;
    real b_dens;
    real b_temp;
    real<lower=0> sigma_proc;
}
transformed parameters{
    vector[n_obs] log_lik;
    vector[n_obs] mu;
    vector[n_day] day_vec;
    for(i in 1:n_day) day_vec[i] = 1;
    for(i in 1:n_obs) {
        mu[i] = b_0[siteID[i]] + b_dens*X_pre[i] + b_temp*W[i,]*day_vec;
        log_lik[i] = normal_lpdf(X_t[i]|mu[i], sigma_proc);
    }
}

model {
    // priors    
    b_0_mean ~ normal(0,2);
    for(j in 1:n_site) {
        b_0[j] ~ normal(b_0_mean, b_0_sigma);
    }
    // lag-1 density coef
    b_dens ~ normal(0,5);
    // temperature coef
    b_temp ~ normal(0,5);
    // residual error
    sigma_proc ~ cauchy(0,5);
    for(i in 1:n_obs) {
            X_t[i] ~ normal(mu[i], sigma_proc);
    }
}
