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
pmcmc = ParticleMarkovChainMonteCarloSV(500, 2500);
pmcmc.run(st);
