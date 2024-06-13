functions {
  vector SpectrumVector(int N,
                        vector lambda,
                        vector lambda_max_vec,
                        real A,
                        real B,
                        real C,
                        real D,
                        real a_1,
                        real a_2,
                        real b,
                        real c,
                        real betaMax_1,
                        real betaMax_2,
                        real bBetaBand_1,
                        real betaBandwidth,
                        real betaBandFactor) {
    vector[N] mu_vec;
    vector[N] lambda_beta_max =  betaMax_1 + betaMax_2 * lambda_max_vec;
    vector[N] b_betaBand = bBetaBand_1 + betaBandwidth .* lambda_max_vec;
    vector[N] beta_band = betaBandFactor * exp( -square((lambda-lambda_beta_max) ./ b_betaBand));
    vector[N] max_div_lambda = lambda_max_vec ./ lambda;
    
    // calculate alpha band with the wavelength (lambda) and the peak (lambda_max_vec) 1/(exp(A(a-x)))
    mu_vec = 1/(exp(A * (a_1+a_2*exp( - square(lambda_max_vec-300)/11940.0) - max_div_lambda))
    + exp(B * (b - max_div_lambda))
    + exp(C * (c - max_div_lambda))
    + D);
    mu_vec .*= exp(A * (a_1+a_2*exp(- square(lambda_max_vec-300)/11940)-1))
    + exp(B * (b-1))
    + exp(C * (c-1))
    + D;
    
    // add beta band and subtract overlap at peak to normalise peak to 1
    mu_vec += beta_band;
    mu_vec ./= 1+betaBandFactor * exp( - square((lambda_max_vec-lambda_beta_max) ./ (b_betaBand)));
    return(mu_vec);
  }
  
  vector LambdaMaxVec(int N,
                      int GroupN,
                      real lambda_max,
                      matrix Group,
                      vector beta_AlphaGroup,
                      matrix betas_Exp,
                      array[] int Experiment) {
    vector[N] lambda_max_vec = exp(log(lambda_max) + Group * beta_AlphaGroup + betas_Exp[Experiment, 3]);
    return(lambda_max_vec);
  }
}

data {
  int<lower=1> N;
  int<lower=0, upper=1> priorOnly; // 1 if only prior is calculated, 0 if likelihood is calculated

  // data to fit
  vector<lower=0>[N] Activity;
  vector<lower=0>[N] lambda;
  vector<lower=0>[N] illumination_time;

  // groups
  int<lower=1> GroupN;
  matrix[N, GroupN-1] Group;
  
  // plates
  int<lower=1> PlateN;
  array[N] int<lower=1> Plate;

  // data for experiment definition
  int<lower=1> ExperimentN;
  array[N] int<lower=1> Experiment;
  array[ExperimentN] int<lower=0> ExperimentConstruct;

  // data for dark activity
  int<lower=1> DarkN;
  int<lower=1> DarkGroupN;
  matrix[DarkN, DarkGroupN-1] Dark_Group;
  array[DarkN] int<lower=1> DarkExperimentID;
  array[DarkN] int<lower=1> DarkPlateID;
  vector<lower=0>[DarkN] DarkActivity;
}

transformed data {
  array[N] int<lower=0> GroupID;
  array[DarkN] int<lower=0> GroupID_dark;
  real lambda_prior_intercept = log(498);
  for(n in 1:N) {
    GroupID[n] = ExperimentConstruct[Experiment[n]];
  }
  
  for(n in 1:DarkN) {
    GroupID_dark[n] = ExperimentConstruct[DarkExperimentID[n]];
  }
}

parameters {
  // error parameter
  real<lower=0> shape;
  
  // template parameters
  // alpha peak parameters
  real<lower=0> A;
  real<lower=0> B;
  real C_raw;
  real c_raw;
  real D_raw;  
  real b_raw;
  real a_1_raw;
  real a_2_raw;
  real lambda_max_raw;
  
  // beta peak parameters
  real<lower=0> betaBandFactor;
  real<lower=0> betaBandwidth;
  real<lower=0> betaMax_1;
  real<lower=0> betaMax_2;
  real bBetaBand_1_raw;
  real<lower=0> betaBandShape;
  real<lower=0> betaBandMean;
  
  // group parameters
  vector[GroupN-1] beta_AlphaGroup_raw;
  vector[GroupN-1] Dark_Group_beta;
  vector[GroupN-1] beta_A_0_Group;
  vector[GroupN-1] beta_A_1_Group;
  
  // dosage-saturation curve parameters
  real<lower=0> A_0;
  real<lower=0> A_1;
  
  // dark activity
  real<lower=0> Dark_shape;
  real Dark_beta;
  
  // hierarchical parameters
  cholesky_factor_corr[3+GroupN] L_Omega_Plate;
  vector<lower=0>[3+GroupN] tau_Plate;  // prior scale
  matrix[3+GroupN, PlateN] z_Plate;
  
  real<lower=0> tau_Group;  // prior scale
  vector[GroupN] z_Group;
}

transformed parameters {
  matrix[PlateN, 3+GroupN] betas_Plate = (diag_pre_multiply(tau_Plate, L_Omega_Plate) * z_Plate)'; // hierarchical coefficients
  vector[GroupN] betas_Group = tau_Group * z_Group; // hierarchical coefficients
  vector[N] Dark_mu_full_log;
  vector[N] Dark_mu_full;

  // construct template parameters estimation
  vector[GroupN-1] beta_AlphaGroup = beta_AlphaGroup_raw * 0.05;
  real lambda_max = exp(lambda_max_raw * 0.07 + lambda_prior_intercept);
  real a_1 = inv_logit(a_1_raw);
  real b = inv_logit(b_raw);
  real a_2 = 0.0459 + a_2_raw * 0.01;
  real C = -14.9 + C_raw * 5;
  real D = 0.674 + D_raw * 0.05;
  real c = 1.104 + c_raw * 0.1; 
  real bBetaBand_1 =  -40.5 + bBetaBand_1_raw * 10;

  // average dark noise for each data point (log scale)
  Dark_mu_full_log = Dark_beta + betas_Plate[Plate, 4] + Group * Dark_Group_beta;
  for(n_add in 2:DarkGroupN) {
    Dark_mu_full_log += Group[,n_add-1] .* betas_Plate[Plate, 3+n_add];
  }
  Dark_mu_full = exp(Dark_mu_full_log);
}

model {
  // template parameters
  target += gamma_lpdf(A | 80, 80/69.7);
  target += gamma_lpdf(B | 20, 20/28.0);
  target += std_normal_lpdf(C_raw);
  target += std_normal_lpdf(D_raw);
  target += normal_lpdf(a_1_raw | 1.987704, 0.45); // mode 0.8795, sd = 0.01
  target += std_normal_lpdf(a_2_raw);
  target += normal_lpdf(b_raw | 2.469836, 0.5); // mode 0.922, sd = 0.041
  target += std_normal_lpdf(c_raw);
  target += gamma_lpdf(betaMax_1 | 500, 500/189.0);
  target += gamma_lpdf(betaMax_2 | 300, 300/0.315);
  target += std_normal_lpdf(bBetaBand_1_raw);
  target += gamma_lpdf(betaBandShape | 5, 0.05);
  target += gamma_lpdf(betaBandMean | 10, 10/0.26);

  target += gamma_lpdf(shape | 5, 5/20.0);
  target += std_normal_lpdf(lambda_max_raw);
  target += gamma_lpdf(betaBandwidth | 20,20/0.195);
  target += gamma_lpdf(betaBandFactor | betaBandShape, betaBandShape/betaBandMean);
  
  // fixed effects
  target += std_normal_lpdf(beta_AlphaGroup_raw);
  target += normal_lpdf(beta_A_0_Group | 0, 0.5);
  target += normal_lpdf(Dark_Group_beta | 0, 1);
  target += normal_lpdf(beta_A_1_Group | 0, .5);

  // dark noise
  target += gamma_lpdf(Dark_shape | 4, 4/20.0);//15.0);
  target += normal_lpdf(Dark_beta | log(0.025), 1);
    
  // dosage-saturation curve
  target += gamma_lpdf(A_0 | 7, 7/0.3);
  target += gamma_lpdf(A_1 | 4, 4/7500.0);
  
  // hierarchical parameters
  target += std_normal_lpdf(to_vector(z_Plate));
  target += lkj_corr_cholesky_lpdf(L_Omega_Plate | 2.0);
  target += normal_lpdf(tau_Plate[1] | 0, 0.1); //shape
  target += normal_lpdf(tau_Plate[2] | 0, 0.05); // A_0
  target += normal_lpdf(tau_Plate[3] | 0, 0.0025); // lambda intercept
  target += normal_lpdf(tau_Plate[4] |0, 0.5); // Dark_mu intercept
  target += normal_lpdf(tau_Plate[5:(3 + GroupN)] | 0, 0.05); // Dark_mu_betas
  
  target += std_normal_lpdf(to_vector(z_Group));
  target += normal_lpdf(tau_Group | 0, 0.5); // shape

  vector[DarkN] Dark_mu;
  vector[N] shape_vec = exp(log(shape) + betas_Plate[Plate, 1] + betas_Group[GroupID]);
  {
    // estimate dark activity
    vector[DarkN] Dark_mu_log = Dark_beta + betas_Plate[DarkPlateID, 4] + Dark_Group * Dark_Group_beta;
    for(n_add in 2:DarkGroupN) {
      Dark_mu_log += Dark_Group[,n_add-1] .* betas_Plate[DarkPlateID, 3+n_add];
    }
    Dark_mu = exp(Dark_mu_log);
  }
  
  vector[N] lambda_max_vec = LambdaMaxVec(N, GroupN, lambda_max, Group, beta_AlphaGroup, betas_Plate, Plate);
  
  // variables for beta band
  vector[N] A1_mu_inv = SpectrumVector(N, lambda, lambda_max_vec, A, B, C, D, a_1, a_2, b, c, betaMax_1, betaMax_2, bBetaBand_1, betaBandwidth, betaBandFactor);
  vector[N] Activity_mu  = exp(log(A_0) + betas_Plate[Plate, 2] + Group * beta_A_0_Group + log(illumination_time)) ./ (exp(log(A_1) + Group * beta_A_1_Group) ./ A1_mu_inv + illumination_time);

  if(!priorOnly) {
    target += gamma_lpdf(DarkActivity | Dark_shape, Dark_shape/Dark_mu);
    target += gamma_lpdf(Activity | shape_vec, shape_vec ./ (Activity_mu+Dark_mu_full));
  }
  
}

generated quantities {
    vector[N] log_lik;
    vector<lower=0>[N] y_rep;
    vector[N] lambda_max_vec = LambdaMaxVec(N, GroupN, lambda_max, Group, beta_AlphaGroup, betas_Plate, Plate);
    vector[N] shape_vec = exp(log(shape) + betas_Plate[Plate, 1] + betas_Group[GroupID]);

    // variables for beta band
    vector[N] A1_mu_inv = SpectrumVector(N, lambda, lambda_max_vec, A, B, C, D, a_1, a_2, b, c, betaMax_1, betaMax_2, bBetaBand_1, betaBandwidth, betaBandFactor);
    vector[N] Activity_mu  = exp(log(A_0) + betas_Plate[Plate, 2] + Group * beta_A_0_Group + log(illumination_time)) ./ (exp(log(A_1) + Group * beta_A_1_Group) ./ A1_mu_inv + illumination_time);

    for(n in 1:N) {
        log_lik[n] = gamma_lpdf(Activity[n] | shape_vec[n], shape_vec[n] / (Activity_mu[n]+Dark_mu_full[n]));
        y_rep[n] = gamma_rng(shape_vec[n],  shape_vec[n]/(Activity_mu[n]+Dark_mu_full[n]));
    }
}
