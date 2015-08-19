function [ vol_two_factors_leverage ] = Sto_Volatility_Estimator( px, steps_mcmc, particles )
    returns = Compute_Returns(px);
    st = SVLogReturns(returns, 0);
    pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(steps_mcmc, particles);
    pmcmc.run(st);
    seq_mcmc = steps_mcmc-500:steps_mcmc;
    rho1    = pmcmc.rho_prop(seq_mcmc);
    sigma1  = pmcmc.sigma_prop(seq_mcmc);
    rho2    = pmcmc.rho2_prop(seq_mcmc);
    sigma2  = pmcmc.sigma2_prop(seq_mcmc);
    beta    = pmcmc.beta_prop(seq_mcmc);
    cor     = pmcmc.cor_prop(seq_mcmc);
    [~, estimated_states_lev, estimated_states2_lev] = BootstrapParticleFilter_TwoFactorsCor(st.y, rho1, sigma1, rho2, sigma2, beta, cor, 10000, pmcmc.p_y_given_x);
    vol_two_factors_leverage = beta^2*exp(estimated_states_lev+estimated_states2_lev);
end

