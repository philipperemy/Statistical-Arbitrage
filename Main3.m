clear;
addpath('filters/'); %rmpath('h/')
addpath('helpers/');
addpath('likelihoods/');
addpath('models/');
addpath('pmcmc/');
addpath('doc/');
diary('log/log.txt');

mySeed = 58; 
%mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

rho = 0.91;       %Scale parameter for X
sigma = 1;        %Standard deviation of the X (process)
beta = 0.5;       %Scale parameter for Y
steps = 1000;     %Number of steps
cor = 0.5;        %Correlation factor between the innovations
mu = 0.2;

st = StochasticVolatilityModelNormalLeverageMean(rho, sigma, beta, steps, cor, mu);
corr(st.innov_X(:), st.innov_Y(:))

%unfortunately. The combination cor = 0.59, mu = 0.60 beats the true values
[estimated_states, return_weights, p_y] = Bootstrap_filter_SIR_NormalLeverageMean(st.y, rho, sigma, beta, 2000, steps, cor, mu);
Y_log_SV = Ex2_Estimate_Py( p_y , steps );
Y_log_SV = Y_log_SV(length(Y_log_SV));
Y_log_SV %-1143.4

%[~, ~, ~, error] = Bootstrap_filter_SIR_NormalLeverage(st.y, -0.593766, 2.705447, beta, 2000, steps, 0.767102)
%Bootstrap_filter_SIR_NormalLeverage(st.y, -0.593766, 2.705447, beta, 2000, steps, 0.767102)
%NaN

pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageMean(2000, 2000);
pmcmc.run(st);
