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

pmcmc = ParticleMarkovChainMonteCarloStudentSV(1, N);
log_val_student = BootstrapParticleFilter_Student(st.y, st.rho, st.sigma, st.beta, st.nu, N, pmcmc.p_y_given_x);

pmcmc = ParticleMarkovChainMonteCarloSV(1, N);
log_val_normal = BootstrapParticleFilter(st.y, st.rho, st.sigma, st.beta, N, pmcmc.p_y_given_x);
bf1 = BayesFactor(log_val_student, log_val_normal);
assert(20 < bf1);

pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverage(1, N);
log_val_normal_leverage_1 = BootstrapParticleFilter_NormalLeverage(st.y, st.rho, st.sigma, st.beta, 0.5, 0, N, pmcmc.p_y_given_x);
log_val_normal_leverage_2 = BootstrapParticleFilter_NormalLeverage(st.y, st.rho, st.sigma, st.beta, 0.7, 0, N, pmcmc.p_y_given_x);

bf2 = BayesFactor(log_val_student, log_val_normal_leverage_1);
bf3 = BayesFactor(log_val_student, log_val_normal_leverage_2);

assert(log_val_normal_leverage_2 < log_val_normal_leverage_1); %2 is less likely than 1
assert(20 < bf2);
assert(20 < bf3);

N = 1000;
pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageMean(1000, N);
pmcmc.run(st);
%mu should be close to 0. cor should also be close to 0. sigma should be
%small
% MCMC_Checks(pmcmc.rho_prop);
% MCMC_Checks(pmcmc.sigma_prop);
% MCMC_Checks(pmcmc.cor_prop);
% MCMC_Checks(pmcmc.mu_prop);
Assert_Range(median(pmcmc.rho_prop), 0.80, 0.95);
Assert_Range(median(pmcmc.sigma_prop), 1.00, 2.00);
Assert_Range(median(pmcmc.cor_prop), 0.00, 0.10);

rho = 0.83;
sigma = 1.52;
cor = 0.08;
mu = 0.60;

N = 100000;
log_val_normal_leverage_mean_1 = BootstrapParticleFilter_NormalLeverage(st.y, 0.83, 1.52, st.beta, 0.08, 0.6, N, pmcmc.p_y_given_x);
bf4 = BayesFactor(log_val_student, log_val_normal_leverage_mean_1);
assert(10 < bf4);