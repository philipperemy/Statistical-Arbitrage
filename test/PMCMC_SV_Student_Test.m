clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');
addpath('../pmcmc/');

load('../data/spx.mat');

mySeed = 57;
rng(mySeed,'twister');

steps_mcmc = 20000;
N = 1000; %1000 is not a lot so good mixing
st = StochasticVolatilityModelStudent;
pmcmc = ParticleMarkovChainMonteCarloStudentSV(steps_mcmc, N);
pmcmc.run(st);
Assert_Range(median(pmcmc.rho_prop), 0.88, 0.93);
Assert_Range(median(pmcmc.sigma_prop), 0.80, 1.20);
Assert_Range(median(pmcmc.nu_prop), 2.50, 3.50);

% N = 3000;
% pmcmc = ParticleMarkovChainMonteCarloStudentSV(st.steps, N);
% pmcmc.run(st);
% MCMC_Checks(pmcmc.rho_prop);
% MCMC_Checks(pmcmc.sigma_prop);
% MCMC_Checks(pmcmc.nu_prop);
% Assert_Range(median(pmcmc.rho_prop), 0.89, 0.92);
% Assert_Range(median(pmcmc.sigma_prop), 0.90, 1.10);
% Assert_Range(median(pmcmc.nu_prop), 2.60, 3.20);


log_val = BootstrapParticleFilter_Student(st.y, st.rho, st.sigma, st.beta, st.nu, N, pmcmc.p_y_given_x);

fprintf('Test completed\n');