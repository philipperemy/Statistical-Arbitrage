classdef ParticleMarkovChainMonteCarloSVNormalLeverage < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop = {};
        sigma_prop = {};
        cor_prop = {};

    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSVNormalLeverage(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.cor_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                %(y, rho_prop, sigma_prop, beta, N, cor)
                case 'rho_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverage(hmm.y, fieldval, this.sigma_prop(step-1), hmm.beta, this.particles, this.cor_prop(step-1));
                case 'sigma_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverage(hmm.y, this.rho_prop(step-1), fieldval, hmm.beta, this.particles, this.cor_prop(step-1));
                case 'cor_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverage(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.particles, fieldval);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = Compute_Log_Likelihood_NormalLeverage(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.particles, this.cor_prop(step-1));
        end
        
        function [val] = call_prior(this, fieldname)
            switch fieldname
                case 'rho_prop'
                    val = -1 + 2*rand; %unif
                case 'sigma_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                case 'cor_prop'
                    val = -1 + 2*rand; %unif
                    %val = 0.9*(-1 + 2*rand); %unif(-0.9, 0.9)
            end
        end 
    end
end