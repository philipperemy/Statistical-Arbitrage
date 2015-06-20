%
% Tests are not very good with this model.
%

clear; clc;
addpath('../filters/');
addpath('../models/');
addpath('../pmcmc/');

mySeed = 57;
rng(mySeed,'twister');

%process 1
rho1 = 0.91;
sigma1 = 1;
%process 2
rho2 = 0.5;
sigma2 = 3;

beta = 0.5;
steps = 1000;
st = StochasticVolatilityModelTwoFactors(rho1, sigma1, rho2, sigma2, beta, steps);

N = 1000;

p_y_given_x = @(y,sig) normpdf(y, 0, sig);
[log_p_y_given_theta, estimated_states, estimated_states2] = BootstrapParticleFilter_TwoFactors(st.y, rho1, sigma1, rho2, sigma2, st.beta, N*100, p_y_given_x);

plot(1:steps, st.x, 'k', 1:steps, estimated_states, 'r');
figure;
plot(1:steps, st.z, 'k', 1:steps, estimated_states2, 'r');

Assert_Range(log_p_y_given_theta, -1691 * 1.01, -1691 * 0.99);

pmcmc = ParticleMarkovChainMonteCarloSVTwoFactors(5000, N*5);
pmcmc.run(st);