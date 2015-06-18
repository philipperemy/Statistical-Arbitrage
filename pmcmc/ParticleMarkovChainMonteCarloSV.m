classdef ParticleMarkovChainMonteCarloSV < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop = {};
        sigma_prop = {};
        beta_prop = {};

        p_y_given_x = @(y,sig) normpdf(y, 0, sig);
    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSV(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.beta_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'rho_prop'
                    log_val = BootstrapParticleFilter(hmm.y, fieldval, this.sigma_prop(step-1), this.beta_prop(step-1), this.particles, this.p_y_given_x);
                case 'sigma_prop'
                    log_val = BootstrapParticleFilter(hmm.y, this.rho_prop(step-1), fieldval, this.beta_prop(step-1), this.particles, this.p_y_given_x);
                case 'beta_prop'
                    log_val = BootstrapParticleFilter(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), fieldval, this.particles, this.p_y_given_x);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.beta_prop(step-1), this.particles, this.p_y_given_x); 
        end
        
        
        function [val] = call_prior(~, fieldname)
            switch fieldname
                case 'rho_prop'
                    val = -1 + 2*rand; %unif
                case 'sigma_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 5)
                       val = 5; 
                    end
                case 'beta_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 5)
                       val = 5; 
                    end
            end
        end
    end
end