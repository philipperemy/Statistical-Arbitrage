clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');
addpath('../coint/deepsearch/');
addpath('../coint/impl/');

%%%% FIRST GET A SPREAD %%%% 09-Sep-2003 - 04-Jun-2006
load('../data/spx.mat');
pp = TruncateData(pp, 5000, 5999);
%spreads = SpreadFinder(pp, 0.8, 236:260);
%AMR CORP	238	 Industrials	 CRANE CO	452	 Industrials	 DOVER CORP	498	 Industrials	0.022396	0.182313	0.776812	0.869268	0.911668	0.833231	0.306703	0.490675	0.353043	1.150421

[ h, spread, ~ ] = SpreadConstructor( [238 452 498], [GetPrice(pp,238) GetPrice(pp,452) GetPrice(pp,498)] );
returns = Compute_Returns(spread.px);
st = SVLogReturns(returns, 0);

steps_mcmc = 5000;
particles = 1000;

subplot(2,1,1);
plot(spread.px, 'b');
xlabel('Time');
ylabel('Price');
subplot(2,1,2);
bar(st.y, 'b');
xlabel('Time');
ylabel('Percentage');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY SIMPLE MODEL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloSV(steps_mcmc, particles);
%pmcmc.run(st);

%MCMC_Checks(pmcmc.rho_prop(1500:5000));
%MCMC_Checks(pmcmc.sigma_prop(1500:5000));
%MCMC_Checks(pmcmc.beta_prop(1500:5000));
%save beta_chain.mat beta_chain
%beta_chain SV std model

%Mean : 0.997985, Median : 0.998154, Min : 0.996784, Max : 0.999437, Acceptance rate ratio : 0.001429 
%Mean : 0.225779, Median : 0.223803, Min : 0.105609, Max : 0.377583, Acceptance rate ratio : 0.114286 
%Mean : 0.697415, Median : 0.441938, Min : 0.125596, Max : 5.000000, Acceptance rate ratio : 0.323714 

%values at the convergence
rho     = 0.998154;
sigma   = 0.223803;
beta 	= 0.441938;

[log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

%yt|xt = N(0, beta^2 exp(xt))
returns_volatility_var = beta^2*exp(estX_SV);
returns_volatility_sd = sqrt(returns_volatility_var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY STUDENT BETA %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloStudentBetaSV(steps_mcmc, particles);
pmcmc.run(st);

MCMC_Checks(pmcmc.beta_prop(2000:5000));
MCMC_Checks(pmcmc.rho_prop(2000:5000));
MCMC_Checks(pmcmc.sigma_prop(2000:5000));
MCMC_Checks(pmcmc.nu_prop(2000:5000));

%Mean : 1.262648, Median : 0.572228, Min : 0.113516, Max : 20.000000, Acceptance rate ratio : 0.455667 
%Mean : 0.991165, Median : 0.999362, Min : -0.959100, Max : 0.999608, Acceptance rate ratio : 0.004000 
%Mean : 0.257189, Median : 0.175191, Min : 0.109505, Max : 19.699607, Acceptance rate ratio : 0.031333 
%Mean : 22.430611, Median : 12.054428, Min : 0.273351, Max : 746.036524, Acceptance rate ratio : 0.113667 

beta = 0.572228;
rho = 0.999362;
sigma = 0.175191;
nu = 12.054428;
[log_p_y_given_theta_SVT, estX_SVT] = BootstrapParticleFilter_Student(st.y, rho, sigma, beta, nu, 10000, pmcmc.p_y_given_x);
BayesFactor(log_p_y_given_theta_SVT, log_p_y_given_theta_SV) %bayes factor is 8-20


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY LEVERAGE BETA %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
particles = 3000;
pmcmc = ParticleMarkovChainMonteCarloSVNormalLeverageBeta(steps_mcmc, particles);
pmcmc.run(st);

%MCMC_Checks(pmcmc.beta_prop(2000:3700));
%MCMC_Checks(pmcmc.rho_prop(2000:3700));
%MCMC_Checks(pmcmc.sigma_prop(2000:3700));
%MCMC_Checks(pmcmc.cor_prop(2000:3700));

%Mean : 0.833002, Median : 0.455972, Min : 0.117021, Max : 20.000000, Acceptance rate ratio : 0.381765 
%Mean : 0.996165, Median : 0.998628, Min : 0.533821, Max : 0.999800, Acceptance rate ratio : 0.005294 
%Mean : 0.239011, Median : 0.218811, Min : 0.118289, Max : 1.500000, Acceptance rate ratio : 0.112353 
%Mean : -0.290313, Median : -0.301758, Min : -0.825591, Max : 0.330991, Acceptance rate ratio : 0.238824 

beta = 0.455972;
rho = 0.998628;
sigma = 0.218811;
cor = -0.301758;
[log_p_y_given_theta_SVL_Beta, estimated_states_SVL] = BootstrapParticleFilter_NormalLeverage(st.y, rho, sigma, beta, cor, 0, 10000, pmcmc.p_y_given_x);
ret_svl_vol = beta^2*exp(estimated_states_SVL)*(1-cor^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY TWO FACTORS     %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = 0.1; %Maybe change beta
st = SVLogReturns(returns, beta);
pmcmc = ParticleMarkovChainMonteCarloSVTwoFactors(steps_mcmc, particles);
pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop(1000:3000));
MCMC_Checks(pmcmc.sigma_prop(1000:3000));
MCMC_Checks(pmcmc.rho2_prop(1000:3000));
MCMC_Checks(pmcmc.sigma2_prop(1000:3000));

%Mean : 0.369315, Median : 0.353178, Min : 0.000756, Max : 0.966032, Acceptance rate ratio : 0.558500 
%Mean : 0.428208, Median : 0.436662, Min : 0.138386, Max : 0.818374, Acceptance rate ratio : 0.142500 
%Mean : 0.997916, Median : 0.997821, Min : 0.996304, Max : 0.998704, Acceptance rate ratio : 0.002000 
%Mean : 0.159545, Median : 0.150783, Min : 0.100266, Max : 0.261483, Acceptance rate ratio : 0.031500 

%MCMC step: 3021 ,rho_prop = 0.110528 ,sigma_prop = 0.453060 ,rho2_prop = 0.999626 ,sigma2_prop = 0.139634 , likelihood = 1802.272854
%stuck infinite loop. M=100

%MCMC step: 1390 ,rho_prop = 0.192676 ,sigma_prop = 0.491330 ,rho2_prop = 0.999559 ,sigma2_prop = 0.126825 ,beta_prop = 0.347858 , likelihood = 1801.176401

rho1 = 0.192676;
sigma1 = 0.491330;
rho2 = 0.999559;
sigma2 = 0.126825;
beta = 0.3478;
[log_p_y_given_theta_TwoF, estimated_states_TwoF, estimated_states_TwoF2] = BootstrapParticleFilter_TwoFactors(st.y, rho1, sigma1, rho2, sigma2, beta, 10000, pmcmc.p_y_given_x);
vol_two_factors = beta^2*exp(estimated_states_TwoF+estimated_states_TwoF2);

Y = sqrt([vol_two_factors' returns_volatility_var' ret_svl_vol']);
plot(Y(10:end,:));
legend('Two Factors SV ', 'Standard SV', 'Normal Leverage SV');

subplot(2,1,1);
plot(estimated_states_TwoF(10:end)');
legend('Latent Process X');
subplot(2,1,2);
plot(estimated_states_TwoF2(10:end)', 'r');
legend('Latent Process Z');

%%% STO TWO FACTORS COR

st = SVLogReturns(returns, 0);
pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(steps_mcmc, particles);
pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop(4500:5000));
MCMC_Checks(pmcmc.sigma_prop(4500:5000));
MCMC_Checks(pmcmc.rho2_prop(4500:5000));
MCMC_Checks(pmcmc.sigma2_prop(4500:5000));
MCMC_Checks(pmcmc.beta_prop(4500:5000));
MCMC_Checks(pmcmc.cor_prop(4500:5000));

%1000-5000
%Mean : 0.726889, Median : 0.939994, Min : 0.000980, Max : 0.999692, Acceptance rate ratio : 0.315000 
%Mean : 0.251178, Median : 0.158994, Min : 0.000055, Max : 1.797519, Acceptance rate ratio : 0.188500 
%Mean : 0.724905, Median : 0.968776, Min : 0.000519, Max : 0.999451, Acceptance rate ratio : 0.273000 
%Mean : 0.293810, Median : 0.223074, Min : 0.009743, Max : 1.256647, Acceptance rate ratio : 0.183250 
%Mean : 0.776786, Median : 0.460051, Min : 0.132099, Max : 5.000000, Acceptance rate ratio : 0.346750 
%Mean : -0.287015, Median : -0.292068, Min : -0.999363, Max : 0.985386, Acceptance rate ratio : 0.341750

%MCMC step: 4950 ,rho_prop = 0.999674 ,sigma_prop = 0.072545 ,rho2_prop = 0.755427 ,sigma2_prop = 0.344350 ,beta_prop = 0.369088 ,cor_prop = -0.853284 , likelihood = 1808.145163

rho1    = 0.999674;
sigma1  = 0.072545;
rho2    = 0.755427;
sigma2  = 0.344350;
beta    = 0.369088;
cor     = -0.853284;

[log_p_y_given_theta_two_factors_lev, estimated_states_lev, estimated_states2_lev] = BootstrapParticleFilter_TwoFactorsCor(st.y, rho1, sigma1, rho2, sigma2, beta, cor, 10000, pmcmc.p_y_given_x);
vol_two_factors_leverage = beta^2*exp(estimated_states_lev+estimated_states2_lev);

%Figure
beg = 10; close all;
figure; plot(estimated_states_lev(beg:end), 'Color', 'blue'); legend('X');
figure;
plot(estimated_states2_lev(beg:end),'Color', 'blue'); legend('Z');
figure;
plot(estimated_states_lev(beg:end)+estimated_states2_lev(beg:end),'Color', 'blue'); legend('X+Z');
figure;
plot(vol_two_factors_leverage(beg:end),'Color', 'blue'); legend('Volatility on returns');
%EndFigure

%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY TWO FACTORS BETA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %beta unknown now OLD SEE WITH beta = 0.3478;
st = SVLogReturns(returns, 0);
pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBeta(steps_mcmc, particles);
pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop(2000:4940));
MCMC_Checks(pmcmc.sigma_prop(2000:4940));
MCMC_Checks(pmcmc.rho2_prop(2000:4940));
MCMC_Checks(pmcmc.sigma2_prop(2000:4940));
MCMC_Checks(pmcmc.beta_prop(2000:4940));

%Mean : 0.370129, Median : 0.367440, Min : 0.000281, Max : 0.966904, Acceptance rate ratio : 0.548639 
%Mean : 0.435403, Median : 0.426310, Min : 0.125813, Max : 0.835795, Acceptance rate ratio : 0.144898 
%Mean : 0.999254, Median : 0.999250, Min : 0.998431, Max : 0.999857, Acceptance rate ratio : 0.001361 
%Mean : 0.164125, Median : 0.164335, Min : 0.099906, Max : 0.311600, Acceptance rate ratio : 0.067347 
%Mean : 1.063325, Median : 0.691514, Min : 0.109365, Max : 5.000000, Acceptance rate ratio : 0.505782 

% rho1 = 0.367440;
% sigma1 = 0.426310;
% rho2 = 0.999250;
% sigma2 = 0.164335;
% beta = 0.691514;
% [log_p_y_given_theta_TwoF, estimated_states_TwoF, estimated_states_TwoF2] = BootstrapParticleFilter_TwoFactors(st.y, rho1, sigma1, rho2, sigma2, beta, 10000, pmcmc.p_y_given_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY MA1             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pmcmc = ParticleMarkovChainMonteCarloSVMA1(steps_mcmc, particles);
pmcmc.run(st);

MCMC_Checks(pmcmc.rho_prop(2000:5000));
MCMC_Checks(pmcmc.sigma_prop(2000:5000));
MCMC_Checks(pmcmc.beta_prop(2000:5000));
MCMC_Checks(pmcmc.psi_prop(2000:5000));

%Mean : 0.999031, Median : 0.999031, Min : 0.999031, Max : 0.999031, Acceptance rate ratio : 0.000000 
%Mean : 0.209434, Median : 0.203622, Min : 0.137682, Max : 0.419117, Acceptance rate ratio : 0.141667 
%Mean : 1.114030, Median : 0.703558, Min : 0.112790, Max : 5.000000, Acceptance rate ratio : 0.533333 
%Mean : 0.483437, Median : 0.476059, Min : 0.000119, Max : 0.999781, Acceptance rate ratio : 0.741667

rho = 0.999031;
sigma = 0.203622;
beta = 0.703558;
psi = 0.476059; %not stable at all
[log_p_y_given_theta_MA1, estimated_states_MA1] = BootstrapParticleFilter_SVMA1(st.y, rho, sigma, beta, psi, 10000, pmcmc.p_y_given_x);

%stuck in an infinite loop... to check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% STOCHASTIC VOLATILITY SVM             %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seems to be very stable. And sigma Acceptance Rate of 15 percent is good.
pmcmc = ParticleMarkovChainMonteCarloSVM(steps_mcmc, particles);
pmcmc.run(st);

%MCMC_Checks(pmcmc.rho_prop(1000:3800));
%MCMC_Checks(pmcmc.sigma_prop(1000:3800));
%MCMC_Checks(pmcmc.beta_prop(1000:3800));

%Mean : 0.998773, Median : 0.998703, Min : 0.997795, Max : 0.999732, Acceptance rate ratio : 0.001786 
%Mean : 0.214640, Median : 0.214255, Min : 0.118362, Max : 0.362811, Acceptance rate ratio : 0.122500 
%Mean : 0.133178, Median : 0.130454, Min : 0.130454, Max : 0.142916, Acceptance rate ratio : 0.000714

rho = 0.998703;
sigma = 0.214255;
beta = 0.130454;
[log_val_SVM, estimated_states_SVM] = BootstrapParticleFilter_SVM(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

save sigma_svm;
load sigma_svm.mat;
sigma_svm = sigma_svm(5000:20000);
sigma_svm(sigma_svm > 0.4) = [];
MCMC_Checks(sigma_svm);

quantile(sigma_svm, 1-0.975);
quantile(sigma_svm, 0.975);
vol_svm = beta^2*exp(estimated_states_SVM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% VOLATILITY MODELLING %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = sqrt([vol_two_factors' returns_volatility_var' ret_svl_vol' vol_two_factors_leverage']);


%%%%%%%%%%% STANDARD SV %%%%%%%%%%%
Volatility_Estimation_Bands_Mean_Zero(spread, returns_volatility_var);
close all;
%%%%%%%%% TWO FACTORS SV (EASIER THAN TFSVL) %%%%%%%%%%%
Volatility_Estimation_Bands_Mean_Zero(spread, vol_two_factors);
close all;
%%%%%%%%% TWO FACTORS SVL %%%%%%%%%%%
Volatility_Estimation_Bands_Mean_Zero(spread, vol_two_factors_leverage);
close all;
Y = [diff_sv diff_tfsv diff_tfsvl];
plot(Y(25:400, :));
legend('SV', 'TFSV', 'TFSVL');