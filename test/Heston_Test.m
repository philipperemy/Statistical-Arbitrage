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

N = 2000;
log_val = BootstrapParticleFilter_Heston(hm.y, mu, theta, kappa, zeta, N);

steps_mcmc = 25000;
N = 2000;
pmcmc = PMCMC_Heston(steps_mcmc, N);
pmcmc.run(hm);



APPLE_TICKER = 239;
appl_prices = GetPrice(pp, APPLE_TICKER);

log_returns = diff(log(appl_prices))';
log_returns = log_returns(5000:5999);

st = SVLogReturns(log_returns, 0);
pmcmc = PMCMC_Heston(steps_mcmc, N);
pmcmc.run(st);

N = 10000;
[log_val, estX] = BootstrapParticleFilter_Heston(st.y, 0.000773, 0.005697, 0.008500, 0.009342, N);

BayesFactor(log_p_y_given_theta_SV, log_val);

%APPL_test
plot([estX_SV(20:end)' estX(20:end)']);
legend('SV Model', 'Heston Model');

