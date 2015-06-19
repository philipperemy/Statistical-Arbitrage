classdef ParticleMarkovChainMonteCarloStudentSV < ParticleMarkovChainMonteCarlo
    
    properties
	
        rho_prop = {};
        sigma_prop = {};
        nu_prop = {};
        
        p_y_given_x = @(y,sig, nu) tpdf(y./sig, nu)./sig;

    end
  
    methods
        function this = ParticleMarkovChainMonteCarloStudentSV(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.rho_prop = zeros(1, this.steps_mcmc);
            this.sigma_prop = zeros(1, this.steps_mcmc);
            this.nu_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'rho_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, fieldval, this.sigma_prop(step-1), hmm.beta, this.nu_prop(step-1), this.particles, this.p_y_given_x);
                case 'sigma_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), fieldval, hmm.beta, this.nu_prop(step-1), this.particles, this.p_y_given_x);
                case 'nu_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, fieldval, this.particles, this.p_y_given_x);
            end
        end
        
        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), hmm.beta, this.nu_prop(step-1), this.particles, this.p_y_given_x); 
        end
        
        
        function [val] = call_prior(~, fieldname)
            switch fieldname
                case 'rho_prop'
                    val = -1 + 2*rand; %unif
                case 'sigma_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                case 'nu_prop'
                    val = 1/gamrnd(1, 1/1); %inverse gamma
            end
        end
        
    end
end