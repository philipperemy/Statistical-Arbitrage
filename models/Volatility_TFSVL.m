function [ vol_two_factors_leverage ] = Volatility_TFSVL( Spread, rho1, sigma1, rho2, sigma2, beta, cor )
    returns = Compute_Returns(Spread.px);
    st = SVLogReturns(returns, 0);
    pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(0, 0);
    [~, estimated_states_lev, estimated_states2_lev] = BootstrapParticleFilter_TwoFactorsCor(st.y, rho1, sigma1, rho2, sigma2, beta, cor, 10000, pmcmc.p_y_given_x);
    vol_two_factors_leverage = beta^2*exp(estimated_states_lev+estimated_states2_lev);
end

