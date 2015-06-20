clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');
addpath('../pmcmc/');

load('../data/spx.mat');

mySeed = 57;
rng(mySeed,'twister');

mu      = 0.01;     %rate of return of the asset 
theta   = 0.02^2;       %long variance
kappa   = 0.2;      %rate at which xt reverts to theta
zeta    = 0.005;        %vol of the vol (variance of xt)
hm = HestonModel(mu, theta, kappa, zeta, 500);

N = 10000;
log_val = BootstrapParticleFilter_Heston(hm.y, mu, theta, kappa, zeta, N);

steps_mcmc = 25000;
N = 10000;
pmcmc = PMCMC_Heston(steps_mcmc, N);
pmcmc.run(hm);
