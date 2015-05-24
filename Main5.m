clear;
addpath('filters/'); %rmpath('h/')
addpath('helpers/');
addpath('likelihoods/');
addpath('models/');
addpath('pmcmc/');
addpath('doc/');
diary('log/log.txt');

st = StochasticVolatilityModelTwoFactors;
pmcmc = ParticleMarkovChainMonteCarloSVTwoFactors(10000, 1000, 'TwoFactors');
% a = pmcmc.run(st);
pmcmc.run(st);