clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');

load('../data/spx.mat');

APPLE_TICKER = 239;
appl_prices = GetPrice(pp, APPLE_TICKER);
appl_prices = appl_prices(5000:5999);
%%% Now we use the returns St/St-1 - 1 instead of the log returns
returns = Compute_Returns(appl_prices);
st = SVLogReturns(returns, 0);

steps_mcmc = 5000;
particles = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY SIMPLE MODEL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
pmcmc = ParticleMarkovChainMonteCarloSV(steps_mcmc, particles);
%pmcmc.run(st);
toc;

%MCMC_Checks(pmcmc.rho_prop(1500:2500));
%MCMC_Checks(pmcmc.sigma_prop(1500:2500));
%MCMC_Checks(pmcmc.beta_prop(1500:2500));

%Mean : 0.998508, Median : 0.999657, Min : 0.424933, Max : 0.999657, Acceptance rate ratio : 0.001001 
%Mean : 0.248592, Median : 0.250599, Min : 0.149552, Max : 5.000000, Acceptance rate ratio : 0.029029 
%Mean : 1.903714, Median : 1.203703, Min : 0.139485, Max : 5.000000, Acceptance rate ratio : 0.690691 
%values at the convergence
rho = 0.9991;
sigma = 0.2353;
beta = 0.8783;

[log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

%yt|xt = N(0, beta^2 exp(xt))
returns_volatility = beta^2*exp(estX_SV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY STUDENT BETA %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloStudentBetaSV(steps_mcmc, particles);
%pmcmc.run(st);

%Mean : 0.994104, Median : 0.998935, Min : -0.970535, Max : 0.998935, Acceptance rate ratio : 0.002000
%Mean : 0.240319, Median : 0.198367, Min : 0.135181, Max : 5.416733, Acceptance rate ratio : 0.016403 
%Mean : 10.154769, Median : 7.685627, Min : 0.609519, Max : 120.077010, Acceptance rate ratio : 0.092218 
%Mean : 0.601529, Median : 0.370564, Min : 0.125216, Max : 14.567718, Acceptance rate ratio : 0.260652 

[log_p_y_given_theta_SVT, estX_SVT] = BootstrapParticleFilter_Student(st.y, 0.9989, 0.1983, 0.3705, 7.685, 10000, pmcmc.p_y_given_x);
BayesFactor(log_p_y_given_theta_SVT, log_p_y_given_theta_SV) %bayes factor is 8-20
% 
% MCMC_Checks(pmcmc.rho_prop(1000:6000));
% figure;
% MCMC_Checks(pmcmc.sigma_prop(1000:6000));
% figure;
% MCMC_Checks(pmcmc.nu_prop(1000:6000));
% figure;
% MCMC_Checks(pmcmc.beta_prop(1000:6000));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY LEVERAGE NO BETA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.1; %add a model to estimate it as well
st = SVLogReturns(returns, beta);
particles = 3000;
pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverage(steps_mcmc, particles);
%pmcmc.run(st);
[log_p_y_given_theta_SVL, estimated_states_SVL] = BootstrapParticleFilter_NormalLeverage(st.y, 0.99586, 0.256, 0.1, -0.435, 0, 10000, pmcmc.p_y_given_x);

ret_svl_vol = beta^2*exp(estimated_states_SVL)*(1-0.435^2);

%MCMC_Checks(pmcmc.rho_prop(1000:5000));
%MCMC_Checks(pmcmc.sigma_prop(1000:5000));
%MCMC_Checks(pmcmc.cor_prop(1000:5000));

%Mean : 0.995290, Median : 0.996000, Min : 0.992173, Max : 0.996000, Acceptance rate ratio : 0.004750 
%Mean : 0.274332, Median : 0.272876, Min : 0.175229, Max : 0.398414, Acceptance rate ratio : 0.099500 
%Mean : -0.433117, Median : -0.439797, Min : -0.622221, Max : -0.123169, Acceptance rate ratio : 0.110000 
%[log_p_y_given_theta_SVL, estimated_states_SVL] = BootstrapParticleFilter_NormalLeverage(st.y, 0.995290, 0.27433, 0.1, -0.4331175, 0, 10000, pmcmc.p_y_given_x);
%Bayes factor is 2.20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY LEVERAGE BETA %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
particles = 3000;
pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageBeta(steps_mcmc, particles);
pmcmc.run(st);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY TWO FACTORS     %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloSVTwoFactors(steps_mcmc, particles);
pmcmc.run(st);

%MCMC step: 807 ,rho_prop = 0.997832 ,sigma_prop = 0.116205 ,rho2_prop = 0.144376 ,sigma2_prop = 0.639814 , likelihood = 2663.636953
%stuck infinite loop. M=50
rho1 = 0.9978;
sigma1 = 0.1162;
rho2 = 0.14437;
sigma2 = 0.6398;
beta = 0.1;
[log_p_y_given_theta_TwoF, estimated_states_TwoF, estimated_states_TwoF2] = BootstrapParticleFilter_TwoFactors(st.y, rho1, sigma1, rho2, sigma2, beta, 10000, pmcmc.p_y_given_x);

vol_two_factors = beta^2*exp(estimated_states_TwoF+estimated_states_TwoF2);

Y = sqrt([vol_two_factors' returns_volatility' ret_svl_vol']);
plot(Y(10:end,:));
legend('Two Factors SV ', 'Standard SV', 'Normal Leverage SV');

subplot(2,1,1);
plot(estimated_states_TwoF(10:end)');
legend('Latent Process X');
subplot(2,1,2);
plot(estimated_states_TwoF2(10:end)', 'r');
legend('Latent Process Z');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY MA1             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloSVMA1(steps_mcmc, particles);
pmcmc.run(st);

%MCMC step: 2749 ,rho_prop = 0.998095 ,sigma_prop = 0.169641 ,beta_prop = 0.235396 ,psi_prop = 0.006020 , likelihood = 2649.250849

rho = 0.998095;
sigma = 0.169641;
beta = 0.235396;
psi = 0.006020;
[log_p_y_given_theta_MA1, estimated_states_MA1] = BootstrapParticleFilter_SVMA1(st.y, rho, sigma, beta, psi, 10000, pmcmc.p_y_given_x);

%stuck in an infinite loop... to check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY SVM             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seems to be very stable. And sigma Acceptance Rate of 15 percent is good.
pmcmc = ParticleMarkovChainMonteCarloSVM(steps_mcmc, particles);
pmcmc.run(st);
MCMC_Checks(pmcmc.rho_prop(1000:5000));
MCMC_Checks(pmcmc.sigma_prop(1000:5000));
MCMC_Checks(pmcmc.beta_prop(1000:5000));
%Mean : 0.998259, Median : 0.998681, Min : 0.996955, Max : 0.998681, Acceptance rate ratio : 0.000500 
%Mean : 0.258502, Median : 0.253392, Min : 0.145597, Max : 0.418673, Acceptance rate ratio : 0.150500 
%Mean : 0.149455, Median : 0.162586, Min : 0.095523, Max : 0.214124, Acceptance rate ratio : 0.003500 

rho = 0.998681;
sigma = 0.253392;
beta = 0.162586;
[log_val_SVM, estimated_states_SVM] = BootstrapParticleFilter_SVM(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

save sigma_svm;
load sigma_svm.mat;
sigma_svm = sigma_svm(5000:20000);
sigma_svm(sigma_svm > 0.4) = [];
MCMC_Checks(sigma_svm);

quantile(sigma_svm, 1-0.975);
quantile(sigma_svm, 0.975);

vol_svm = exp(estimated_states_SVM);

Y = sqrt([vol_two_factors' returns_volatility' ret_svl_vol' vol_svm']);
plot(Y(10:end,:));