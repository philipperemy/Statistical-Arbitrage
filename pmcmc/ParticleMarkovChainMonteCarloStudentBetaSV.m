classdef ParticleMarkovChainMonteCarloStudentBetaSV < ParticleMarkovChainMonteCarloStudentSV
    
    properties
        %Check in base class
        beta_prop = {};
    end
  
    methods
        function this = ParticleMarkovChainMonteCarloStudentBetaSV(steps_mcmc, particles)
            this = this@ParticleMarkovChainMonteCarloStudentSV(steps_mcmc, particles);
            this.beta_prop = zeros(1, this.steps_mcmc);
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'beta_prop' 
                    log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), fieldval, this.nu_prop(step-1), this.particles, this.p_y_given_x);
                case 'rho_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, fieldval, this.sigma_prop(step-1), this.beta_prop(step-1), this.nu_prop(step-1), this.particles, this.p_y_given_x);
                case 'sigma_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), fieldval, this.beta_prop(step-1), this.nu_prop(step-1), this.particles, this.p_y_given_x);
                case 'nu_prop'
                    log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.beta_prop(step-1), fieldval, this.particles, this.p_y_given_x);
            end
        end

        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_Student(hmm.y, this.rho_prop(step-1), this.sigma_prop(step-1), this.beta_prop(step-1), this.nu_prop(step-1), this.particles, this.p_y_given_x); 
        end
        
        function [val] = call_prior(this, fieldname)
            switch fieldname
                case 'beta_prop' 
                    val = 1/gamrnd(1, 1/1); %inverse gamma
                    if(val > 20)
                       val = 20; 
                    end
                otherwise
                    val = call_prior@ParticleMarkovChainMonteCarloStudentSV(this, fieldname);
            end
        end
        
    end
end