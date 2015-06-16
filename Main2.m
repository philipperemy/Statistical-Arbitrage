clear;
addpath('filters/');
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

st = StochasticVolatilityModelNormalLeverage(rho, sigma, beta, steps, cor);
corr(st.innov_X(:), st.innov_Y(:))
% 
% seq = -0.90:0.05:0.90;
% correlations = size(seq,2);
% i = 1;
% for c = seq
%     [estimated_states, return_weights, p_y] = Bootstrap_filter_SIR_NormalLeverage(st.y, rho, sigma, beta, 2000, steps, c);
%     Y_log_SV = Ex2_Estimate_Py( p_y , steps );
%     Y_log_SV = Y_log_SV(length(Y_log_SV));
%     correlations(i) = Y_log_SV;
%     i = i + 1;
% end
% plot(seq, correlations);
% seq(correlations == max(correlations))
% 
% plot(1:steps, st.x, 'k', 1:steps, estimated_states, 'r');

[estimated_states, return_weights, p_y] = Bootstrap_filter_SIR_NormalLeverage(st.y, rho, sigma, beta, 2000, steps, cor);
Y_log_SV = Ex2_Estimate_Py( p_y , steps );
Y_log_SV = Y_log_SV(length(Y_log_SV));
Y_log_SV %-1063.9

%[~, ~, ~, error] = Bootstrap_filter_SIR_NormalLeverage(st.y, -0.593766, 2.705447, beta, 2000, steps, 0.767102)
%Bootstrap_filter_SIR_NormalLeverage(st.y, -0.593766, 2.705447, beta, 2000, steps, 0.767102)
%NaN

pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverage(2000, 2000);
pmcmc.run(st);
