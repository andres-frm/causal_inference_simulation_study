functions {
  matrix GP_quadratic(matrix x, // distance matrix
                      real eta, 
                      real rho, 
                      real delta) {
                        
                        int N = dims(x)[1];
                        matrix[N, N] K;
                        
                        for (i in 1:(N-1)) {
                          K[i, i] = eta + delta;
                          for (j in (i+1):N) {
                            K[i, j] = eta * exp(-rho*square(x[i, j]));
                            K[j, i] = K[i, j];
                          }
                        }
                        K[N,N] = eta + delta;
                        return(cholesky_decompose(K));
                      }
}

data {
  int N;
  int N_study;
  int N_sites;
  int N_spp;
  int N_genus;
  int N_family;
  int N_na;
  int N_no_na;
  array[N] int study;
  array[N] int sites;
  array[N] int sp;
  array[N] int genus;
  array[N] int family;
  array[N_na] int index_na;
  array[N_no_na] int index_no_na;
  vector[N] X;
  vector[N] Y;
  vector[N_no_na] fruit_availability_OBS;
  matrix[N_study, N_study] dist_mat;
  
}

parameters {
  
  // random effects
  vector[N_study] z_study_Y;
  real<lower = 0> rho_Y;
  real<lower = 0> eta_Y;
  
  vector[N_sites] z_sites_Y;
  real mu_sites_Y;
  real<lower = 0> sigma_sites_Y;
  
  vector[N_study] z_study_X;
  real<lower = 0> rho_X;
  real<lower = 0> eta_X;
  
  vector[N_study] z_study_f;
  real<lower = 0> rho_f;
  real<lower = 0> eta_f;
  
  vector[N_sites] z_sites_X;
  real mu_sites_X;
  real<lower = 0> sigma_sites_X;
  
  vector[N_sites] z_sites_F;
  real mu_sites_F;
  real<lower = 0> sigma_sites_F;
  
  vector[N_spp] z_sp;
  real mu_sp;
  real<lower = 0> sigma_sp;
  
  vector[N_genus] z_genus;
  real mu_genus;
  real<lower = 0> sigma_genus;
  
  vector[N_family] z_family;
  real mu_family;
  real<lower = 0> sigma_family;
  
  vector[N_spp] z_sp_x;
  real mu_sp_x;
  real<lower = 0> sigma_sp_x;
  
  vector[N_genus] z_genus_x;
  real mu_genus_x;
  real<lower = 0> sigma_genus_x;
  
  vector[N_family] z_family_x;
  real mu_family_x;
  real<lower = 0> sigma_family_x;
  
  vector[N_na] fruit_missing;
  
  real beta_fruit_X;
  real beta_fruit_Y;
  real beta_X;
  real alpha_X;
  real alpha_f;
  real alpha_Y;
  real<lower = 0> sigma_Y;
  real<lower = 0> sigma_f;
  real<lower = 0> sigma_X;
  
}

transformed parameters{
  matrix[N_study, N_study] K_Ly;
  vector[N_study] STUDY_Y;
  K_Ly = GP_quadratic(dist_mat, eta_Y, rho_Y, 0.001);
  STUDY_Y = K_Ly * z_study_Y;
  
  matrix[N_study, N_study] K_Lx;
  vector[N_study] STUDY_X;
  K_Lx = GP_quadratic(dist_mat, eta_X, rho_X, 0.001);
  STUDY_X = K_Lx * z_study_X;
  
  matrix[N_study, N_study] K_Lf;
  vector[N_study] STUDY_F;
  K_Lf = GP_quadratic(dist_mat, eta_f, rho_f, 0.001);
  STUDY_F = K_Lf * z_study_f;
  
  vector[N_sites] SITE_Y;
  SITE_Y = mu_sites_Y + z_sites_Y * sigma_sites_Y;
  
  vector[N_sites] SITE_X;
  SITE_X = mu_sites_X + z_sites_X * sigma_sites_X;
  
  vector[N_sites] SITE_F;
  SITE_F = mu_sites_F + z_sites_F * sigma_sites_F;
  
  vector[N_spp] SP;
  SP = mu_sp + z_sp * sigma_sp;
  
  vector[N_genus] GENUS;
  GENUS = mu_genus + z_genus * sigma_genus;
  
  vector[N_family] FAMILY;
  FAMILY = mu_family + z_family * sigma_family;
  
  vector[N_spp] SP_x;
  SP_x = mu_sp_x + z_sp_x * sigma_sp_x;
  
  vector[N_genus] GENUS_x;
  GENUS_x = mu_genus_x + z_genus_x * sigma_genus_x;
  
  vector[N_family] FAMILY_x;
  FAMILY_x = mu_family_x + z_family_x * sigma_family_x;
  
  vector[N] fruit_all;
  fruit_all[index_no_na] = fruit_availability_OBS;
  fruit_all[index_na] = fruit_missing;
  
}

model {
  alpha_f ~ normal(0, 1); 
  alpha_X ~ normal(0, 1); 
  alpha_Y ~ normal(0, 1); 
  
  beta_X ~ normal(0, 1);
  
  beta_fruit_X ~ normal(0, 1);
  
  beta_fruit_Y ~ normal(0 , 1);
  
  eta_Y ~ exponential(3);
  z_study_Y ~ normal(0, 0.5);
  rho_Y ~ exponential(1);
  
  eta_X ~ exponential(3);
  z_study_X ~ normal(0, 0.5);
  rho_X ~ exponential(1);
  
  eta_f ~ exponential(3);
  z_study_f ~ normal(0, 0.5);
  rho_f ~ exponential(1);
  
  z_sites_Y ~ normal(0, 0.5);
  mu_sites_Y ~ normal(0, 0.25);
  sigma_sites_Y ~ exponential(2);
  
  z_sites_X ~ normal(0, 0.5);
  mu_sites_X ~ normal(0, 0.25);
  sigma_sites_X ~ exponential(2);
  
  z_sites_F ~ normal(0, 0.5);
  mu_sites_F ~ normal(0, 0.25);
  sigma_sites_F ~ exponential(2);
  
  z_sp ~ normal(0, 0.5);
  mu_sp ~ normal(0, 0.25);
  sigma_sp ~ exponential(2);
  
  z_genus ~ normal(0, 0.5);
  mu_genus ~ normal(0, 0.25);
  sigma_genus ~ exponential(2);
  
  z_family ~ normal(0, 1);
  mu_family ~ normal(0, 0.25);
  sigma_family ~ exponential(2);
  
  z_sp_x ~ normal(0, 1);
  mu_sp_x ~ normal(0, 0.25);
  sigma_sp_x ~ exponential(2);
  
  z_genus_x ~ normal(0, 1);
  mu_genus_x ~ normal(0, 0.25);
  sigma_genus_x ~ exponential(2);
  
  z_family_x ~ normal(0, 0.5);
  mu_family_x ~ normal(0, 0.25);
  sigma_family_x ~ exponential(2);
  
  sigma_f ~ exponential(1);
  sigma_Y ~ exponential(1);
  sigma_X ~ exponential(1);
  
  fruit_all ~ normal(alpha_f + STUDY_F[study] + SITE_F[sites],
                     sigma_f);
  
  X ~ normal(alpha_X + STUDY_X[study] + SITE_X[sites] +
                  SP_x[sp] + GENUS_x[genus] + FAMILY_x[family] +
                  beta_fruit_X*fruit_all, 
                  sigma_X);
  
  Y ~ normal(alpha_Y + STUDY_Y[study] + SITE_Y[sites] +
                  SP[sp] + GENUS[genus] + FAMILY[family] +
                  beta_X*X + beta_fruit_Y*fruit_all, 
                  sigma_Y);
  
}

generated quantities {
  array[N] real ppcheck;

  ppcheck = normal_rng(alpha_Y + STUDY_Y[study] + SITE_Y[sites] +
                       SP[sp] + GENUS[genus] + FAMILY[family] +
                       beta_X*X + beta_fruit_Y*fruit_all,
                       sigma_Y);
}
