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
appl_prices = appl_prices(5000:5999);
%%% Now we use the returns St/St-1 - 1 instead of the log returns
returns = Compute_Returns(appl_prices);
st = SVLogReturns(returns, 0);
steps_mcmc = 5000;
particles = 1000;

%%%%%%%%%%%% STOCHASTIC VOLATILITY SIMPLE MODEL %%%%%%%%%%%%
tic;
pmcmc = ParticleMarkovChainMonteCarloSV(steps_mcmc, particles);
%pmcmc.run(st);
toc;

% MCMC_Checks(pmcmc.rho_prop(1:1000));
% MCMC_Checks(pmcmc.sigma_prop(1:1000));
% MCMC_Checks(pmcmc.beta_prop(1:1000));

%Mean : 0.998508, Median : 0.999657, Min : 0.424933, Max : 0.999657, Acceptance rate ratio : 0.001001 
%Mean : 0.248592, Median : 0.250599, Min : 0.149552, Max : 5.000000, Acceptance rate ratio : 0.029029 
%Mean : 1.903714, Median : 1.203703, Min : 0.139485, Max : 5.000000, Acceptance rate ratio : 0.690691 
%values at the convergence
rho = 0.9985;
sigma = 0.2505;
beta = 1.203;

[log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

%yt|xt = N(0, beta^2 exp(xt))
returns_volatility_var = beta^2*exp(estX_SV);
returns_volatility_sd = sqrt(returns_volatility_var);
N = 1000; %generate N observations of the process
T = length(returns);
generated_processes = zeros(N, T);

real_vol = zeros(T,1);
for t = 2:T
    real_vol(t) = returns_volatility_sd(t)*appl_prices(t-1);
    generated_processes(:,t) = normrnd(appl_prices(t-1), returns_volatility_sd(t)*appl_prices(t-1), N, 1);
    %generated_processes(:,t) = appl_prices(t-1)*(1+returns_volatility_sd(t)*randn(N,1));
end

plot(generated_processes');

ma = tsmovavg(appl_prices,'s',20,1);

boll_plus = ma + 2*real_vol;
boll_minus = ma - 2*real_vol;
sum((boll_minus < appl_prices) & (appl_prices < boll_plus))
Y = [boll_plus boll_minus appl_prices];

std_1 = movingstd(appl_prices, 3, 'backward');
plot([real_vol std_1]);
%as period->0, we find the same volatility measure

%do it just with two MA
boll_plus_2 = zeros(T,1);
boll_minus_2 = zeros(T,1);
ma = tsmovavg(appl_prices,'s',5,1);
real_vol_ma = tsmovavg(real_vol,'s',5,1);
for t = 2:T
    boll_plus_2(t) = ma(t-1) + 1.5*real_vol_ma(t);
    boll_minus_2(t) = ma(t-1) - 1.5*real_vol_ma(t);
end
Y = [boll_plus_2 boll_minus_2 appl_prices];


%%%No moving average. Only prices. Very sensitive.
boll_plus_3 = zeros(T,1);
boll_minus_3 = zeros(T,1);
for t = 2:T
    boll_plus_3(t) = appl_prices(t-1) + real_vol(t);
    boll_minus_3(t) = appl_prices(t-1) - real_vol(t);
end
Y = [boll_plus_3 boll_minus_3 appl_prices];

%%%%Second method - use generated process
movstd_mat = zeros(N,T);
for n = 1:N
    movstd_mat(n,:) = movingstd(generated_processes(n,:), 20, 'backward');
end

plot(movstd_mat');
figure;
plot(movingstd(appl_prices, 20, 'backward'));
plot([mean(movstd_mat, 1)' movingstd(appl_prices, 20, 'backward')]);
plot([quantile(movstd_mat, 0.9, 1)' movingstd(appl_prices, 20, 'backward')]);

plot(mean(movstd_mat, 1));
hold on;
plot(movingstd(appl_prices, 20, 'backward'));

std_eval = mean(movstd_mat, 1)';
std_real = movingstd(appl_prices, 20, 'backward');
ma = tsmovavg(appl_prices,'s',20,1);
d = 1.5;
Y = [ma-d*std_eval ma appl_prices ma+d*std_eval];
Y = [Y ma-d*std_real ma+d*std_real];
plot(Y);
legend('boll-e', 'ma', 'prices', 'boll+e', 'boll-r', 'boll+r');

normalized_std_diff = (std_eval-std_real)./appl_prices;
plot(normalized_std_diff(30:end))