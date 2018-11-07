data {
    int N; // number of observations
    int N2; // number of generated observations
    vector[N] q; // observations
    vector[N] years_obs; // years of observations
    vector<lower=0>[N] sigma; // standard deviations
    vector[N2] years; // years to predict
}
parameters {
  real <lower=1970> t_m; // time when mdr arose
  real <lower=0> b; // slope
  real rho; // dummy parameter
}
transformed parameters {
  real c;
  vector[N] quadpred;
  vector[N] yrs_fromtm;

  c = rho * b / t_m;
  yrs_fromtm = years_obs - rep_vector(t_m,N);
  
  for(i in 1:N)
    quadpred[i] = b * yrs_fromtm[i] - c * pow(yrs_fromtm[i],2);
}
model { 
  q ~ normal(quadpred, sigma);
  t_m ~ normal(1985, 9);
  b ~ lognormal(-5.5, 0.7);
  rho ~ normal(5,15)T[,36]; 
}
generated quantities {
  vector[N2] p_pred;
  vector[N2] yrsp_fromtm;
  
  yrsp_fromtm = years - rep_vector(t_m,N2);
  
  for(i in 1:N2)
    p_pred[i] = b * yrsp_fromtm[i] - c * pow(yrsp_fromtm[i],2); //observations predicted by the model
}
