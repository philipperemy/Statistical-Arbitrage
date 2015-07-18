classdef ParticleMarkovChainMonteCarloSVNormalLeverage < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop    = {};
        sigma_prop  = {};
        cor_prop    = {};
        
        p_y_given_x = @(y,sig) normpdf(y, 0, sig);

    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSVNormalLeverage(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop   = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.cor_prop   = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'rho_prop'
                    log_val = BootstrapParticleFilter_NormalLeverage(hmm.y, fieldval, this.sigma_prop(step-1), hmm.beta, this.cor_prop(step-1), 0, this.particles, this.p_y_given_x);
                case 'sigma_prop'
                    log_val = BootstrapParticleFilter_NormalLeverage(hmm.y, this.rho_prop(step-1), fieldval, hmm.beta, this.cor_prop(step-1), 0, this.particles, this.p_y_given_x);
                case 'cor_prop'
                    log_val = BootstrapParticleFilter_NormalLeverage(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, fieldval, 0, this.particles, this.p_y_given_x);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_NormalLeverage(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.cor_prop(step-1), 0, this.particles, this.p_y_given_x);
        end
        
        function [val] = call_prior(~, fieldname)
            switch fieldname
                case 'rho_prop'
                    val = rand; %unif
                    if(val > 0.996)  %it's because it crashes when -> 1
                        val = 0.996;
                    end
                case 'sigma_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 1.5)
                       val = 1.5; 
                    end
                case 'cor_prop'
                    val = -1 + 2*rand; %unif
                    %val = rand;
            end
        end 
    end
end