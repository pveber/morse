functions {
  real[] TKTD( real t, // time
               real[] y, // variable
               real[] param, // parameters
               real[] x_r, // real data
               int[] x_i) { // integer data

    // - parameters
    int n_group;

    real kd;
    real kk;
    real hb;
    real z;

    // - new variables
    real max_z[2 ,x_i[2]]; // [2, double nb group]
    
    real dy_dt[x_i[2]]; // DOUBLE NUMBER OF GROUP

    n_group = x_i[1];

    // - parameters allocation
    hb = param[1];
    kd = param[2];
    kk = param[3];
    z =  param[4];

    // - model
    for( i in 1:n_group){
      dy_dt[i] = kd * (x_r[i] - y[i]);
    }

    for( i in 1:n_group){
      max_z[1,i] = 0;
      max_z[2,i] = y[i] - z;
    }

    for( i in 1:n_group){
       dy_dt[i+n_group]  = kk * max(max_z[,i]) + hb;
    }

    return dy_dt;
  }

}

data {

  int<lower=1> Tm;
  int<lower=1> n_group;

  real t0;
  real ts[Tm];

  // OBSERVATION
  real conc_cst[n_group];

  int Nsurv[Tm,n_group];
  int Nprec[Tm,n_group];

  // PRIORS
  real kk_meanlog10;
  real kk_sdlog10;
  real z_meanlog10;
  real z_sdlog10;
  real kd_meanlog10;
  real kd_sdlog10;
  real hb_meanlog10;
  real hb_sdlog10;

}
transformed data {
  
  real x_r[n_group];
  int x_i[2];
  
  int twice_n_group;
  
  twice_n_group = 2*n_group;
  x_i[1] = n_group;
  x_i[2] = twice_n_group;

  x_r = conc_cst;
  //x_r = conc_cst;

}
parameters {
  real<lower=0> y0[twice_n_group]; // Must be positive, and length is DOUBLE of n_group

  real kk_log10;
  real z_log10;
  real kd_log10;
  real hb_log10;

}

transformed parameters{

  real<lower=0> param[4]; //

  param[1] = 10^hb_log10; // hb
  param[2] = 10^kd_log10; // kd
  param[3] = 10^kk_log10; // kk
  param[4] = 10^z_log10; // z

}
model {

  real y_hat[Tm,twice_n_group]; // DOUBLE nb of group [Tm, 2*n_group]
  real psurv[Tm,n_group]; // [Tm, n_group]

   kk_log10 ~ normal( kk_meanlog10, kk_sdlog10 );
   z_log10  ~ normal( z_meanlog10,   z_sdlog10 );
   kd_log10 ~ normal( kd_meanlog10, kd_sdlog10 );
   hb_log10 ~ normal( hb_meanlog10, hb_sdlog10 );

  //y0 ~ exponential(10^6); // Initial condition for y0 have to be put close to 0 !!!

  y_hat = integrate_ode_rk45(TKTD, y0, t0, ts, param, x_r, x_i);

   for(j in 1:n_group){
    psurv[1,j] = exp( -y_hat[1, j+n_group]);
    Nsurv[1,j] ~ binomial( Nprec[1,j], psurv[1,j] / 1);

    for(t in 2:Tm){
     psurv[t,j] = exp( -y_hat[t,j+n_group]);
     Nsurv[t,j] ~ binomial( Nprec[t,j], psurv[t,j] / psurv[t-1,j]);
    }
  }
}
generated quantities {
  real y_tilde[Tm,twice_n_group]; // DOUBLE nb of group [Tm, 2*n_group]
  real psurv_tilde[Tm,n_group]; // [Tm, n_group]
  int Nsurv_tilde[Tm,n_group]; // [Tm, n_group]

  y_tilde = integrate_ode_rk45(TKTD, y0, t0, ts, param, x_r, x_i);

  for(j in 1:n_group){
    psurv_tilde[1,j] = exp( -y_tilde[1, j+n_group]);
    Nsurv_tilde[1,j] = binomial_rng( Nprec[1,j], psurv_tilde[1,j] / 1);

    for(t in 2:Tm){
     psurv_tilde[t,j] = exp( -y_tilde[t,j+n_group]);
     Nsurv_tilde[t,j] = binomial_rng( Nsurv_tilde[t-1,j], psurv_tilde[t,j] / psurv_tilde[t-1,j]);
    }
  }
}

