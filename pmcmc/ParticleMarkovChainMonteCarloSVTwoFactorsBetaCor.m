classdef ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop    = {};
        sigma_prop  = {};
        rho2_prop   = {};
        sigma2_prop = {};
        beta_prop   = {};
		cor_prop	= {};

        p_y_given_x = @(y,sig) normpdf(y, 0, sig);
    end
  
    methods
        function this = ParticleMarkovChainMonteCarloSVTwoFactorsBetaCor(steps_mcmc, particles)
            this             = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop    = zeros(1, this.steps_mcmc);
            this.sigma_prop  = zeros(1, this.steps_mcmc);
            this.rho2_prop   = zeros(1, this.steps_mcmc);
            this.sigma2_prop = zeros(1, this.steps_mcmc);
            this.beta_prop   = zeros(1, this.steps_mcmc);
			this.cor_prop	 = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'rho_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, fieldval, this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), this.beta_prop(step-1), this.cor_prop(step-1), this.particles, this.p_y_given_x);
                case 'sigma_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), fieldval, this.rho2_prop(step-1), this.sigma2_prop(step-1), this.beta_prop(step-1), this.cor_prop(step-1), this.particles, this.p_y_given_x);
                case 'rho2_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), fieldval, this.sigma2_prop(step-1), this.beta_prop(step-1), this.cor_prop(step-1), this.particles, this.p_y_given_x);
                case 'sigma2_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), fieldval, this.beta_prop(step-1), this.cor_prop(step-1), this.particles, this.p_y_given_x);
                case 'beta_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), fieldval, this.cor_prop(step-1), this.particles, this.p_y_given_x);
				case 'cor_prop'
                    log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), this.beta_prop(step-1), fieldval, this.particles, this.p_y_given_x);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_TwoFactorsCor(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.rho2_prop(step-1), this.sigma2_prop(step-1), this.beta_prop(step-1), this.cor_prop(step-1), this.particles, this.p_y_given_x); 
        end
        
        function [val] = call_prior(~, fieldname)
            switch fieldname
                case 'beta_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 5)
                       val = 5; 
                    end
                case {'rho_prop'}
                   val = 0.90 + (1-0.90)*rand;
                    % val = rand; %unif, must be positive
                case {'rho2_prop'}
                    val = rand;
                case {'sigma_prop', 'sigma2_prop'}
                    val = 2*rand;
				case 'cor_prop'
					val = 2*rand-1;
                    %val = -rand;
            end
        end
        
    end
end