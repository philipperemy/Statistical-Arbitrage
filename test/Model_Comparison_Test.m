clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');
addpath('../pmcmc/');
addpath('../likelihoods/');

load('../data/spx.mat');

mySeed = 57;
rng(mySeed,'twister');

N = 100000;
st = StochasticVolatilityModelStudent;
pmcmc = ParticleMarkovChainMonteCarloStudentSV(1, N); %just for pmcmc.p_y_given_x
log_val_student = BootstrapParticleFilter_Student(st.y, st.rho, st.sigma, st.beta, st.nu, N, pmcmc.p_y_given_x);

pmcmc = ParticleMarkovChainMonteCarloSV(1, N); %just for pmcmc.p_y_given_x
log_val_normal = BootstrapParticleFilter(st.y, st.rho, st.sigma, st.beta, N, pmcmc.p_y_given_x);
bf = BayesFactor(log_val_student, log_val_normal);
assert(bf > 20);

% st = normal
% estimate it with a t-student dist. nu should be very high.
st = StochasticVolatilityModel;
N = 3000;
pmcmc = ParticleMarkovChainMonteCarloStudentSV(1000, N);
pmcmc.run(st);
MCMC_Checks(pmcmc.rho_prop);
MCMC_Checks(pmcmc.sigma_prop);
MCMC_Checks(pmcmc.nu_prop);

Assert_Range(mean(pmcmc.rho_prop), 0.89, 0.92);
Assert_Range(mean(pmcmc.sigma_prop), 0.95, 1.05);
Assert_Range(mean(pmcmc.nu_prop), 20, 1000);

% st = still normal
% Bayes Factor should tell: normal model is better than t-student model for
% low nu
N = 100000;

pmcmc_student = ParticleMarkovChainMonteCarloStudentSV(1, N); %just for pmcmc.p_y_given_x
pmcmc_normal = ParticleMarkovChainMonteCarloSV(1, N); %just for pmcmc.p_y_given_x

nu = 3;
log_val_student = BootstrapParticleFilter_Student(st.y, st.rho, st.sigma, st.beta, nu, N, pmcmc_student.p_y_given_x);
log_val_normal = BootstrapParticleFilter(st.y, st.rho, st.sigma, st.beta, N, pmcmc_normal.p_y_given_x);
assert(log_val_student < log_val_normal);
bf = BayesFactor(log_val_normal, log_val_student);
assert(bf > 20);

nu = 1000;
log_val_student_1000 = BootstrapParticleFilter_Student(st.y, st.rho, st.sigma, st.beta, nu, N, pmcmc_student.p_y_given_x);
bf = BayesFactor(log_val_normal, log_val_student_1000);
assert(bf < 2); %models are equivalent
aic_student = 2*3 - 2*log_val_student_1000;
aic_normal = 2*2 - 2*log_val_normal;
assert(aic_normal < aic_student);