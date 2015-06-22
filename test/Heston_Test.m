clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');
addpath('../pmcmc/');
addpath('../likelihoods/');

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
%pmcmc.run(hm);



APPLE_TICKER = 239;
appl_prices = GetPrice(pp, APPLE_TICKER);

log_returns = diff(log(appl_prices))';
log_returns = log_returns(5000:5999);

st = SVLogReturns(log_returns, 0);
pmcmc = PMCMC_Heston(steps_mcmc, N);
%pmcmc.run(st);

N = 10000;
[log_val, estX_Heston] = BootstrapParticleFilter_Heston(st.y, 0.000773, 0.005697, 0.008500, 0.009342, N);

BayesFactor(log_p_y_given_theta_SV, log_val);

%APPL_test
%same heston and SV if

% [log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);
% estX_SV = exp(estX_SV); %ok. See. Should have beta^2

plot([estX_SV(20:end)' estX_Heston(20:end)']);
legend('SV Model', 'Heston Model');

a10 = exp(movingstd(log_returns, 10));
a5 = movingstd(log_returns, 5);

a5_b = exp(movingstd(log_returns, 5,'b'));
a5_f = exp(movingstd(log_returns, 5,'f'));

%std and not variance!
b = beta*exp(estX_SV/2); %ok. See. Should have beta^2

%b = exp(estX_SV); %lost! maybe forever. shit. Will find it again for sure.
plot([b(20:end)' a5(20:end)']);
legend('SV', 'std');


log_prices = log(appl_prices);
log_prices = log_prices(5000:5999); %not exact. but still good to have an idea
ma = tsmovavg(log_prices', 's', 20);

b_plus = ma'+200*b.^2';
b_minus = ma'-200*b.^2';
plot([b_plus b_minus log_prices ])

[mid, low, up] = bollinger(log_prices);
plot([low up log_prices ])