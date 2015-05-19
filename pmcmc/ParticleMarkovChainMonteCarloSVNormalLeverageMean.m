classdef ParticleMarkovChainMonteCarloSVNormalLeverageMean < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop = {};
        sigma_prop = {};
        cor_prop = {};
        mu_prop = {};

    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSVNormalLeverageMean(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.cor_prop = zeros(1, this.steps_mcmc);
            this.mu_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                %Compute_Log_Likelihood_NormalLeverageMean(y, rho_prop, sigma_prop, beta, N, cor, mu)
                case 'rho_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverageMean(hmm.y, fieldval, this.sigma_prop(step-1), hmm.beta, this.particles, this.cor_prop(step-1), this.mu_prop(step-1));
                case 'sigma_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverageMean(hmm.y, this.rho_prop(step-1), fieldval, hmm.beta, this.particles, this.cor_prop(step-1), this.mu_prop(step-1));
                case 'cor_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverageMean(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.particles, fieldval, this.mu_prop(step-1));
                case 'mu_prop'
                    log_val = Compute_Log_Likelihood_NormalLeverageMean(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.particles, this.cor_prop(step-1), fieldval);
            
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = Compute_Log_Likelihood_NormalLeverageMean(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.particles, this.cor_prop(step-1), this.mu_prop(step-1));
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
                case 'mu_prop'
                    val = -1 + 2*rand; %unif
            end
        end 
    end
end