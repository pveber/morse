functions {
  /* bisectioning search of the index i such that sorted[i] < x < sorted[i+1]
  this come from Sebastian Weber
  */
    int find_interval_elem(real x, vector sorted, int start_ind) {
      int res;
      int N;
      int max_iter;
      real left;
      real right;
      int left_ind;
      int right_ind;
      int iter;

      N = num_elements(sorted);

      if(N == 0) return(0);

      left_ind  = start_ind;
      right_ind = N;

      max_iter = 100 * N;
      left  = sorted[left_ind ] - x;
      right = sorted[right_ind] - x;

      if(0 <= left)  return(left_ind-1);
      if(0 == right) return(N-1);
      if(0 >  right) return(N);

      iter = 1;
      while((right_ind - left_ind) > 1  && iter != max_iter) {
        int mid_ind;
        real mid;
        // is there a controlled way without being yelled at with a warning?
        mid_ind = (left_ind + right_ind) / 2;
        mid = sorted[mid_ind] - x;
        if (mid == 0) return(mid_ind-1);
        if (left  * mid < 0) { right = mid; right_ind = mid_ind; }
        if (right * mid < 0) { left  = mid; left_ind  = mid_ind; }
        iter = iter + 1;
      }
      if(iter == max_iter)
        print("Maximum number of iterations reached.");
      return(left_ind);
    }

  real linearInterp( real t_x, real t_before, real t_after, real y_before, real y_after){
    real linInterp_hat;

    linInterp_hat = y_before + (t_x - t_before) * (y_after - y_before)/(t_after-t_before);

    return(linInterp_hat);
  }

  real[] varSD_MULTIgroup( real t,
                           real[] y,
                           real[] param,
                           real[] x_r,
                           int[] x_i) {

    // - parameters
    int n_group = x_i[1];
    int Nconc[n_group] = x_i[2:n_group+1];
    int Nconc_cumsum[n_group] = x_i[n_group+2:2*n_group+2];

    // - parameters allocation
    real kd = param[1];


    real dy_dt[n_group]; // DOUBLE NUMBER OF GROUP

    // for variable concentration

    for(gr in 1:n_group){

      vector[Nconc[gr]] tconc_edo = to_vector(x_r[Nconc_cumsum[gr]:Nconc_cumsum[gr+1]]); // to take dose time which is in ts
      int d = find_interval_elem(t, tconc_edo, 1);

      vector[Nconc[gr]] conc_edo = to_vector(x_r[Nconc_cumsum[gr]+1:2*Nconc_cumsum[gr+1]]);

      real linInterp_conc;

      linInterp_conc = linearInterp(t , tconc_edo[d], tconc_edo[d+1], conc_edo[d], conc_edo[d+1] );

      // - model
      dy_dt[gr] =  kd * ( linInterp_conc - y[gr]);
    }

    return dy_dt;
  }
}
data {
  int<lower=1> Tm;
  int<lower=1> n_group; // e.g. 4 = N_profile

  int<lower=1> N_data_conc; // e.g. 6+6+6+6 sum(Nconc)
  int<lower=1> Nconc[n_group]; // e.g. 6 6 6 6
  int<lower=1> Nconc_cumsum[n_group+1]; // e.g. 6 12 18 24
    
  real t0;
  real ts[Tm];

  real conc[N_data_conc]; // size = max(Nconc) ; e.g. 24
  real tconc[N_data_conc]; // size = max(Nconc) ; e.g. 24

  // OBSERVATION
  int Nsurv[Tm,n_group];
  int Nprec[Tm,n_group];

  // PRIORS
  real kd_meanlog10;
  real kd_sdlog10;
  real hb_meanlog10;
  real hb_sdlog10;

  real alpha_meanlog10;
  real alpha_sdlog10;
  real beta_minlog10;
  real beta_maxlog10;

}
transformed data {
  int x_i[2*n_group+2];
  real x_r[2*N_data_conc];

  x_i[1] = n_group;
  x_i[2:n_group+1] = Nconc;
  x_i[n_group+2:2*n_group+2] = Nconc_cumsum;

  x_r[1:N_data_conc] = tconc; // all concentration
  x_r[N_data_conc+1:2*N_data_conc] = conc; // all time of concentration

}
parameters {
  real<lower=0> y0[n_group]; // Must be positive, and length is DOUBLE of n_group

  real kd_log10;
  real hb_log10;
  real alpha_log10;
  real beta_log10;

}

transformed parameters{

  real<lower=0> param[1]; //

  real<lower=0> hb;
  real<lower=0> alpha;
  real beta;

  param[1] = 10^kd_log10; // kd

  hb = 10^hb_log10;
  alpha = 10^alpha_log10;
  beta = 10^beta_log10;
}
model {
  real y_hat[Tm,2*n_group]; // DOUBLE nb of group [Tm, 2*n_group]

  real D_max[Tm,n_group];
  real F_prob[Tm,n_group];
  real psurv_hat[Tm,n_group]; // [Tm, n_group]

  kd_log10 ~ normal( kd_meanlog10, kd_sdlog10 );
  hb_log10 ~ normal( hb_meanlog10, hb_sdlog10 );

  alpha_log10  ~ normal( alpha_meanlog10, alpha_sdlog10 );
  beta_log10 ~ uniform( beta_minlog10 , beta_maxlog10 );

  //y0 ~ exponential(10^6); // Initial condition for y0 have to be put close to 0 !!!

    y_hat = integrate_ode_rk45(varSD_MULTIgroup, y0, t0, ts, param, x_r, x_i);

  for(j in 1:n_group){
    # LOG LOGISTIC
    D_max[1,j] = y_hat[1,j]; // Peut-être trouver le max entre ts=0 et ts=1 ?
    F_prob[1,j] = D_max[1,j]^beta/(D_max[1,j]^beta+alpha^beta);

    psurv_hat[1,j] =  exp(-hb * ts[1]) * (1- F_prob[1,j]);
    Nsurv[1,j] ~ binomial( Nprec[1,j], psurv_hat[1,j]/1);

    for(t in 2:Tm){
      D_max[t,j] = max(y_hat[1:t,j]);
      psurv_hat[t,j] = exp(-hb * ts[t]) * (1- F_prob[t,j]);
      Nsurv[t,j] ~ binomial( Nprec[t,j], psurv_hat[t,j] / psurv_hat[t-1,j]);
    }
  }

}
