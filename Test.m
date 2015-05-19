mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

diary('log.txt');

clear;

st = StochasticVolatilityModel;

pmcmc = ParticleMarkovChainMonteCarloSV(500, 2500);
% a = pmcmc.run(st);
pmcmc.run(st);
