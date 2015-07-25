clear;
addpath('../filters/');
addpath('../models/');
addpath('../helpers/');
addpath('../pmcmc/');

mySeed = 57;
rng(mySeed,'twister');

st = StochasticVolatilityModelM;
pmcmc = ParticleMarkovChainMonteCarloSVM(2000, 1000);
pmcmc.run(st);

Assert_Range(median(pmcmc.rho_prop), 0.88, 0.93); 
Assert_Range(median(pmcmc.sigma_prop), 0.98, 1.02);
Assert_Range(median(pmcmc.beta_prop), 0.40, 0.60);

%parallel may be a little worst for rho 0.884 instead of 0.910
%so it's better to remove it