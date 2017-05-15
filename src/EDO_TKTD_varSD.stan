functions {
  real[] TKTD( real t, // time
               real[] y, // variable
               real[] prm, // parameters
               real[] x_r, // real data
               int[] x_i) { // integer data
    
    // - parameters  
    real ke;
    real kk;
    real hb;
    real z;
    
    real[]  C;
    
    // - new variables
    real max_z[2]; // [2, double nb group]
    real dydt[2]; // DOUBLE NUMBER OF GROUP

    
    // - parameters allocation
    hb = prm[1];
    ke = prm[2];
    kk = prm[3];
    z =  prm[4];
    

    // - linear interpolation   
    C = C[ts[i]] + (t + ts[i]) * (C(ts[i+1]) - C(ts[i]))/(ts[i+1] - ts[i]);
    
    
    dydt[1] = ke * (C - y[1]);
    
    max_z[1] = 0;
    max_z[2] = y[1] - z;
    
    dydt[2]  = kk * max(max_z) + hb;

    return dydt;
  }
  
  
  // real[] lin_interpol(real[] x,
  //                     real[] x0,
  //                     real[] x1,
  //                     real[] y0,
  //                     real[] y1){
  //   
  //   real[]  y;
  //   
  //   y = y0 + (x + x0) * (y1 - y0)/(x1 - x0);
  //   
  //   return y;
  // }
    
}

data {
  
  int<lower=1> Tm;
  int<lower=1> n_group;
  
  real t0;
  real ts[Tm];
  
  // OBSERVATION
  real Conc[Tm,n_group];
  //real conc_cst[n_group];
  
  int Nsurv[Tm,n_group];
  int Nprec[Tm,n_group];
  
  // PRIORS
  real meanlog10_kk;
  real taulog10_kk;
  real meanlog10_z;
  real taulog10_z;
  real meanlog10_ke;
  real taulog10_ke;
  real meanlog10_hb;
  real taulog10_hb;
  
}
transformed data {
  int twice_n_group;
  
  real x_r[n_group];
  int x_i[Tm];
  

  //x_r = Conc[1,];
  x_r = Conc[,1];
  
}
parameters {
  real<lower=0> y0[12]; // Must be positive, and length is DOUBLE of n_group
  
  real log10_kk; 
  real log10_z; 
  real log10_ke; 
  real log10_hb; 
  
}

transformed parameters{
  
  real<lower=0> prm[4]; // kk, z, ke, hb, n_group
  
  prm[1] = 10^log10_hb; // hb
  prm[2] = 10^log10_ke; // ke
  prm[3] = 10^log10_kk; // kk
  prm[4] = 10^log10_z; // z
  
}

model {

  real y_hat[Tm,12]; // DOUBLE nb of group [Tm, 2*n_group]
  real psurv_hat[Tm,6]; // [Tm, n_group]
  
  log10_kk ~ normal(meanlog10_kk, sqrt(1/taulog10_kk));
  log10_z ~ normal(meanlog10_z, sqrt(1/taulog10_z));
  log10_ke ~ normal(meanlog10_ke, sqrt(1/taulog10_ke));
  log10_hb ~ normal(meanlog10_hb, sqrt(1/taulog10_hb));
  
  y0 ~ exponential(10^6); // Initial condition for y0 have to be put close to 0 !!!
  
  y_hat = integrate_ode_rk45(TKTD, y0, t0, ts, prm, x_r, x_i);
  
  for(t in 1:Tm){
    
    for(j in 1:n_group){
       psurv_hat[t,j] = exp(- y_hat[t,j+n_group]);
       Nsurv[t,j] ~ binomial(Nprec[t,j], psurv_hat[t,j]);
    }
  }
  
}