clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');

load('../data/spx.mat');

p_y_given_x = @(y,sig) normpdf(y, 0, sig);

APPLE_TICKER = 239;
appl_prices = GetPrice(pp, APPLE_TICKER);

log_returns = diff(log(appl_prices))';
log_returns = log_returns(5000:5999);

st = SVLogReturns(log_returns, 0);

steps_mcmc = 50000;
particles = 1000;

%%%%%%%%%%%% STOCHASTIC VOLATILITY SIMPLE MODEL %%%%%%%%%%%%
tic;
pmcmc = ParticleMarkovChainMonteCarloSV(20, particles);
pmcmc.run(st);
toc;

rho = 0.9989;
sigma = 0.2243;
beta = 0.7502;

[log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);
estX_SV = exp(estX_SV);
%%%%%%%%%%%% STOCHASTIC VOLATILITY STUDENT MODEL %%%%%%%%%%%%
%assuming beta. Will relax after
beta = 0.7502;
st = SVLogReturns(log_returns, beta);
pmcmc = ParticleMarkovChainMonteCarloStudentSV(steps_mcmc, particles);
%pmcmc.run(st);

%MCMC_Checks(pmcmc.rho_prop(1:1000));
%MCMC_Checks(pmcmc.sigma_prop(1:1000));
%MCMC_Checks(pmcmc.nu_prop(1:1000));

rho     =   0.9989;
sigma   =   0.240;
nu      =   9;

[log_p_y_given_theta_SVT, ~] = BootstrapParticleFilter_Student(st.y, rho, sigma, beta, nu, 10000, pmcmc.p_y_given_x);

BayesFactor(log_p_y_given_theta_SVT, log_p_y_given_theta_SV)

%Relaxing beta condition. estimate it now.
steps_mcmc = 50000;
pmcmc = ParticleMarkovChainMonteCarloStudentBetaSV(steps_mcmc, particles*3);
pmcmc.run(st);

[log_p_y_given_theta_SVT, ~] = BootstrapParticleFilter_Student(st.y, 0.9989, 0.266967, 0.998863, 6.74, 10000, pmcmc.p_y_given_x);
BayesFactor(log_p_y_given_theta_SVT, log_p_y_given_theta_SV) %bayes factor is 8-20

MCMC_Checks(pmcmc.rho_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.sigma_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.nu_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.beta_prop(1000:6000));
