clear;
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
pmcmc.run(st);

tic;
[Y_log_SV, estX] = Ex4_Compute_Log_Likelihood(st.y, st.rho, st.sigma, st.beta, 10000, st.nu);
toc;
[Y_log_SV2, estX2] = Ex4_Compute_Log_Likelihood(st.y, st.rho, st.sigma, st.beta, 10000, st.nu);
plot([estX' estX2' st.x']);
legend('ESS','Simple Res', 'orig');

plot([estX' estX2']);
legend('ESS','Simple Res');

%
% No ESS
%MCMC_Checks(pmcmc.rho_prop(1:8600))
%Mean : 0.878422, Median : 0.903982, Min : -0.987159, Max : 0.992850, Acceptance rate ratio : 0.041633 
%MCMC_Checks(pmcmc.sigma_prop(1:8600))
%Mean : 5.821610, Median : 1.108391, Min : 0.194667, Max : 9440.838710, Acceptance rate ratio : 0.091755 
%MCMC_Checks(pmcmc.nu_prop(1:8600))
%Mean : 5.681863, Median : 2.780274, Min : 0.213094, Max : 2407.856001, Acceptance rate ratio : 0.143156 

%
% ESS stratified
%MCMC_Checks(pmcmc.rho_prop(1:8600))
%Mean : 0.861349, Median : 0.897350, Min : -0.972093, Max : 0.995569, Acceptance rate ratio : 0.049075 
%MCMC_Checks(pmcmc.sigma_prop(1:8600))
%Mean : 33.802678, Median : 1.149064, Min : 0.140347, Max : 213439.214158, Acceptance rate ratio : 0.093267 
%MCMC_Checks(pmcmc.nu_prop(1:8600))
%Mean : 113.785169, Median : 2.852150, Min : 0.301428, Max : 44854.253305, Acceptance rate ratio : 0.144319

%
% ESS multinomial

