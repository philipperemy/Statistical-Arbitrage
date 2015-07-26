clear;
addpath('../filters/');
addpath('../models/');
addpath('../helpers/');
addpath('../pmcmc/');

mySeed = 57;
rng(mySeed,'twister');

st = StochasticVolatilityModelMA1;
pmcmc = ParticleMarkovChainMonteCarloSVMA1(2000, 1000);
pmcmc.run(st);

Assert_Range(median(pmcmc.rho_prop), 0.88, 0.93); 
Assert_Range(median(pmcmc.sigma_prop), 0.90, 1.02);
Assert_Range(median(pmcmc.beta_prop), 0.45, 0.55);
Assert_Range(median(pmcmc.psi_prop), 0.75, 0.85);

log_val = BootstrapParticleFilter_SVMA1(st.y, st.rho, st.sigma, st.beta, st.psi, 10000, pmcmc.p_y_given_x);