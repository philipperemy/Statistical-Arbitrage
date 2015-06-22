classdef PMCMC_Heston2 < ParticleMarkovChainMonteCarlo
    
    properties
        alpha_prop      = {};          %rate of return of the asset 
        beta_prop       = {};          %long variance
        sigma_prop      = {};          %rate at which xt reverts to theta
        theta_prop      = {};          %vol of the vol (variance of xt)
        r_prop          = {};
    end
  
    methods
        function this = PMCMC_Heston2(steps_mcmc, particles)
            this               = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.alpha_prop    = zeros(1, this.steps_mcmc);
            this.beta_prop     = zeros(1, this.steps_mcmc);
            this.sigma_prop    = zeros(1, this.steps_mcmc);
            this.theta_prop    = zeros(1, this.steps_mcmc);
            this.r_prop        = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            alpha = this.alpha_prop(step-1);
            beta = this.beta_prop(step-1);
            sigma = this.sigma_prop(step-1);
            theta = this.theta_prop(step-1);
            r = this.r_prop(step-1);
            
            switch fieldname  
                case 'alpha_prop'
                    alpha = fieldval; 
                case 'beta_prop'
                    beta = fieldval;
                case 'sigma_prop'
                    sigma = fieldval;
                case 'theta_prop'
                    theta = fieldval;
                case 'r_prop'
                    r = fieldval;
            end
            
            log_val = BootstrapParticleFilter_Heston2(hmm.y, alpha, beta, sigma, theta, r, this.particles);  
        end

        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_Heston2(hmm.y, this.alpha_prop(step-1), this.beta_prop(step-1), this.sigma_prop(step-1), this.theta_prop(step-1), this.r_prop(step-1), this.particles); 
        end

        
        function [val] = call_prior(~, ~)
            val = rand;
        end
        
    end
end