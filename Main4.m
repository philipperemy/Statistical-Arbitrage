clear;
addpath('filters/'); %rmpath('h/')
addpath('helpers/');
addpath('likelihoods/');
addpath('models/');
addpath('pmcmc/');
addpath('doc/');
diary('log/log.txt');

mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

rho = 0.91;       %Scale parameter for X
sigma = 1;        %Standard deviation of the X (process)

rho2 = 0.5;
sigma2 = 3;

beta = 0.5;       %Scale parameter for Y

st = StochasticVolatilityModelTwoFactors;
steps = st.steps;     %Number of steps

particles = 3500;

[Y_log_SV] = Compute_Log_Likelihood_TwoFactors(st.y, rho, sigma, rho2, sigma2, beta, particles);
Y_log_SV

[estimated_states, estimated_states2, return_weights, p_y] = Bootstrap_filter_SIR_TwoFactors(st.y, rho, sigma, rho2, sigma2, beta, particles, steps);

plot(1:steps, st.x, 'k', 1:steps, estimated_states, 'r');
plot(1:steps, st.z, 'k', 1:steps, estimated_states2, 'r');

pmcmc = ParticleMarkovChainMonteCarloSVTwoFactors(10000, particles);
% a = pmcmc.run(st);
pmcmc.run(st);
