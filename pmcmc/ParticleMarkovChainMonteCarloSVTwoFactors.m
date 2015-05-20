classdef ParticleMarkovChainMonteCarloSVTwoFactors < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop = {};
        sigma_prop = {};
        
        rho2_prop = {};
        sigma2_prop = {};

    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSVTwoFactors(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.rho2_prop = zeros(1, this.steps_mcmc);
            this.sigma2_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'rho_prop'
                    log_val = Compute_Log_Likelihood_TwoFactors(hmm.y, fieldval, this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), hmm.beta, this.particles);
                case 'sigma_prop'
                    log_val = Compute_Log_Likelihood_TwoFactors(hmm.y, this.rho_prop(step-1), fieldval, this.rho2_prop(step-1), this.sigma2_prop(step-1), hmm.beta, this.particles);
                case 'rho2_prop'
                    log_val = Compute_Log_Likelihood_TwoFactors(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), fieldval, this.sigma2_prop(step-1), hmm.beta, this.particles);
                case 'sigma2_prop'
                    log_val = Compute_Log_Likelihood_TwoFactors(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), fieldval, hmm.beta, this.particles);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = Compute_Log_Likelihood_TwoFactors(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), hmm.beta, this.particles); 
        end
        
        function [val] = call_prior(this, fieldname)
            switch fieldname
                case {'rho_prop', 'rho2_prop'}
                    val = -1 + 2*rand; %unif
                case {'sigma_prop', 'sigma2_prop'}
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 10) %little proba of having two consecutive vals > 10
                        val = 1/gamrnd(1, 1/1); %inverse gamma
                    end
            end
        end
        
    end
end