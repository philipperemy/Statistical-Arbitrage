clear; clc;
addpath('filters/');
addpath('helpers/');
addpath('likelihoods/');
addpath('models/');
addpath('pmcmc/');
addpath('doc/');
diary('log.txt');

mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

st = StochasticVolatilityModelStudent;
pmcmc = ParticleMarkovChainMonteCarloSV(25000, 1000);

RHO = 'rho_prop';
SIG = 'sigma_prop';
NU = 'nu_prop';
MAX_STEP = 1000;

N = 1000;

rho_p_best = 0;
sig_p_best = 10;
nu_p_best  = 10;
denom = Ex4_Compute_Log_Likelihood_ESS(st.y, rho_p_best, sig_p_best, st.beta, N, nu_p_best);

%Monolithic update
for step = 2:MAX_STEP
    
    rho_p = rho_p_best + normrnd(0, 0.2);
    num = Ex4_Compute_Log_Likelihood_ESS(st.y, rho_p, sig_p_best, st.beta, N, nu_p_best);
    if(num > denom)
        rho_p_best = rho_p;
        denom = num;
    end
    clear rho_p;
    
    sig_p = sig_p_best + normrnd(0, 0.2);
    num = Ex4_Compute_Log_Likelihood_ESS(st.y, rho_p_best, sig_p, st.beta, N, nu_p_best);
    if(num > denom)
        sig_p_best = sig_p;
        denom = num;
    end
    clear sig_p;
    
    nu_p = nu_p_best + normrnd(0, 0.2);
    num = Ex4_Compute_Log_Likelihood_ESS(st.y, rho_p_best, sig_p_best, st.beta, N, nu_p);
    if(num > denom)
        nu_p_best = nu_p;
        denom = num;
    end
    clear nu_p;

    fprintf('[INFO] MCMC step: %i ,rho_prop = %f, sigma_prop = %f, nu_prop = %f , likelihood = %f\n',...
         step, rho_p_best, sig_p_best, nu_p_best, num);
    
end

%RHO
% rho_proposal    = normrnd(pmcmc.rho_prop(step-1), 0.1);%pmcmc.call_prior(RHO);
% sig_proposal    = normrnd(pmcmc.sigma_prop(step-1), 0.1);%pmcmc.call_prior(SIG);
% nu_proposal     = 2+normrnd(pmcmc.nu_prop(step-1), 0.1);%pmcmc.call_prior(NU)
% 
% denom = pmcmc.call_likelihood_denom(st, step);
% num = Ex4_Compute_Log_Likelihood_ESS(st.y, rho_proposal, sig_proposal, st.beta, N, nu_proposal);
% 
% fprintf('[DEBUG] num = %f, denom = %f, rho_p = %f, sig_p = %f, nu_p = %f\n', num, denom, rho_proposal, sig_proposal, nu_proposal);
% 
% if(num > denom)
%     fprintf('[INFO] MCMC step: %i ,rho_prop = %f, sigma_prop = %f, nu_prop = %f , likelihood = %f\n',...
%         step, rho_proposal, sig_proposal, nu_proposal, num);
%     pmcmc.rho_prop(step) = rho_proposal;
%     pmcmc.sigma_prop(step) = sig_proposal;
%     pmcmc.nu_prop(step) = nu_proposal;
% else
%     pmcmc.rho_prop(step) = pmcmc.rho_prop(step-1);
%     pmcmc.sigma_prop(step) = pmcmc.sigma_prop(step-1);
%     pmcmc.nu_prop(step) = pmcmc.nu_prop(step-1);
% end