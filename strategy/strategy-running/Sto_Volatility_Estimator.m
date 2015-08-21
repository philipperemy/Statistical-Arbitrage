function [ res ] = Sto_Volatility_Estimator( px, steps_mcmc, particles )
    returns = Compute_Returns(px);
    st = SVLogReturns(returns, 0);
    
    while(true)
        pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(steps_mcmc, particles);
        pmcmc.run(st);
        if(steps_mcmc < 500)
            seq_mcmc = 1:steps_mcmc;
        else
            seq_mcmc = steps_mcmc-500:steps_mcmc;
        end
        rho1    = median(pmcmc.rho_prop(seq_mcmc));
        sigma1  = median(pmcmc.sigma_prop(seq_mcmc));
        rho2    = median(pmcmc.rho2_prop(seq_mcmc));
        sigma2  = median(pmcmc.sigma2_prop(seq_mcmc));
        beta    = median(pmcmc.beta_prop(seq_mcmc));
        cor     = median(pmcmc.cor_prop(seq_mcmc));

        if(rho1 ~= 0)
            break;
        end
    end
    
    fprintf('rho1 = %f, sigma1 = %f, rho2 = %f, sigma2 = %f, beta = %f, cor = %f\n', rho1, sigma1, rho2, sigma2, beta, cor);
    [~, estimated_states_lev, estimated_states2_lev] = BootstrapParticleFilter_TwoFactorsCor(st.y, rho1, sigma1, rho2, sigma2, beta, cor, 10000, pmcmc.p_y_given_x);
    vol_two_factors_leverage = beta^2*exp(estimated_states_lev+estimated_states2_lev);
    res = struct('vol_two_factors_leverage', vol_two_factors_leverage, 'pmcmc', pmcmc );
end

