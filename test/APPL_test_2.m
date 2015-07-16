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

%%% Now we use the returns St/St-1 - 1 instead of the log returns
returns = Compute_Returns(appl_prices);
returns = returns(5000:5999);

st = SVLogReturns(returns, 0);

steps_mcmc = 5000;
particles = 1000;

%%%%%%%%%%%% STOCHASTIC VOLATILITY SIMPLE MODEL %%%%%%%%%%%%
tic;
pmcmc = ParticleMarkovChainMonteCarloSV(steps_mcmc, particles);
pmcmc.run(st);
toc;

MCMC_Checks(pmcmc.rho_prop(1:1000));
MCMC_Checks(pmcmc.sigma_prop(1:1000));
MCMC_Checks(pmcmc.beta_prop(1:1000));

%Mean : 0.998508, Median : 0.999657, Min : 0.424933, Max : 0.999657, Acceptance rate ratio : 0.001001 
%Mean : 0.248592, Median : 0.250599, Min : 0.149552, Max : 5.000000, Acceptance rate ratio : 0.029029 
%Mean : 1.903714, Median : 1.203703, Min : 0.139485, Max : 5.000000, Acceptance rate ratio : 0.690691 
%values at the convergence
rho = 0.9985;
sigma = 0.2505;
beta = 1.203;

[log_p_y_given_theta_SV, estX_SV] = BootstrapParticleFilter(st.y, rho, sigma, beta, 10000, pmcmc.p_y_given_x);

%yt|xt = N(0, beta^2 exp(xt))
returns_volatility = beta^2*exp(estX_SV);


%Relaxing beta condition. estimate it now.
steps_mcmc = 5000;
pmcmc = ParticleMarkovChainMonteCarloStudentBetaSV(steps_mcmc, particles);
pmcmc.run(st);

[log_p_y_given_theta_SVT, estX_SVT] = BootstrapParticleFilter_Student(st.y, 0.9989, 0.266967, 0.998863, 6.74, 10000, pmcmc.p_y_given_x);
BayesFactor(log_p_y_given_theta_SVT, log_p_y_given_theta_SV) %bayes factor is 8-20

MCMC_Checks(pmcmc.rho_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.sigma_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.nu_prop(1000:6000));
figure;
MCMC_Checks(pmcmc.beta_prop(1000:6000));

%define beta before! does not work for student T! but let's do it anyway
returns_volatility_student = beta^2*exp(estX_SVT);
st.x = returns_volatility_student;
N = 1000; %generate N observations of the process
T = length(st.y);
generated_processes = zeros(N, T);

for t = 2:T
    generated_processes(:,t) = normrnd(st.y(t-1), returns_volatility_student(t), N, 1);
end

