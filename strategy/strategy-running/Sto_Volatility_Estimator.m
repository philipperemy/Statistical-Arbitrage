function [ res ] = Sto_Volatility_Estimator( px, steps_mcmc )
    returns = Compute_Returns(px);
    st = SVLogReturns(returns, 0);
    
    particles = length(returns)+100;
    
    i = 1;
    while(true)
        pmcmc = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(steps_mcmc, particles);
        pmcmc.run(st);
        rho1    = pmcmc.rho_prop(end);
        sigma1  = pmcmc.sigma_prop(end);
        rho2    = pmcmc.rho2_prop(end);
        sigma2  = pmcmc.sigma2_prop(end);
        beta    = pmcmc.beta_prop(end);
        cor     = pmcmc.cor_prop(end);

        if(rho1 ~= 0)
            break;
        end
        
        if(i > 2)
           res = struct('vol_two_factors_leverage', 0, 'pmcmc', 0 );
           return;
        end
        
        i = i + 1;
    end
    
    fprintf('rho1 = %f, sigma1 = %f, rho2 = %f, sigma2 = %f, beta = %f, cor = %f\n', rho1, sigma1, rho2, sigma2, beta, cor);
    [~, estimated_states_lev, estimated_states2_lev] = BootstrapParticleFilter_TwoFactorsCor(st.y, rho1, sigma1, rho2, sigma2, beta, cor, 5000, pmcmc.p_y_given_x);
    vol_two_factors_leverage = beta^2*exp(estimated_states_lev+estimated_states2_lev);
    
    pmcmc = struct('rho1', rho1, 'sigma1', sigma1, 'rho2', rho2, 'sigma2', sigma2, 'beta', beta, 'cor', cor);
    res = struct('vol_two_factors_leverage', vol_two_factors_leverage, 'pmcmc', pmcmc );
end

