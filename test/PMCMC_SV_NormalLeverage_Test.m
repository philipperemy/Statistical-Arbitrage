clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');
addpath('../pmcmc/');

load('../data/spx.mat');

mySeed = 57;
rng(mySeed,'twister');

% mu = 0
N = 1000; %1000 is not a lot so good mixing
st = StochasticVolatilityModelNormalLeverage;
Assert_Range(corr(st.innov_Y', st.innov_X'), 0.67, 0.73); %cor = 0.7 (default)

pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverage(st.steps, N);
[true_val, estX] = BootstrapParticleFilter_NormalLeverage(st.y, st.rho, st.sigma, st.beta, st.cor, 0, 10000, pmcmc.p_y_given_x);

pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop);
MCMC_Checks(pmcmc.sigma_prop);
MCMC_Checks(pmcmc.cor_prop);

Assert_Range(median(pmcmc.rho_prop), 0.89, 0.92);
Assert_Range(median(pmcmc.sigma_prop), 0.90, 1.10);
Assert_Range(median(pmcmc.cor_prop), 0.5, 0.8);

% mu = 0.2
st = StochasticVolatilityModelNormalLeverage(st.rho, st.sigma, st.beta, st.steps, st.cor, 0.2);
Assert_Range(corr(st.innov_Y', st.innov_X'), 0.67, 0.73); %cor = 0.7 (default)

pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageMean(st.steps, N);
pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop);
MCMC_Checks(pmcmc.sigma_prop);
MCMC_Checks(pmcmc.cor_prop);
MCMC_Checks(pmcmc.mu_prop);

Assert_Range(median(pmcmc.rho_prop), 0.89, 0.92);
Assert_Range(median(pmcmc.sigma_prop), 0.90, 1.10);
Assert_Range(median(pmcmc.cor_prop), 0.55, 0.85);
Assert_Range(median(pmcmc.mu_prop), -0.5, 0.5); % acc rate = 0.343343 big
%so maybe not precise enough

%N = 3000;
%pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageMean(st.steps, N);
%pmcmc.run(st);
%[true_val, estX] = BootstrapParticleFilter_NormalLeverage(st.y, st.rho, st.sigma, st.beta, st.cor, st.mu, 10000, pmcmc.p_y_given_x);