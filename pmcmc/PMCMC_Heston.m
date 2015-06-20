classdef PMCMC_Heston < ParticleMarkovChainMonteCarlo
    
    properties
        mu_prop      = {};          %rate of return of the asset 
        theta_prop   = {};          %long variance
        kappa_prop   = {};          %rate at which xt reverts to theta
        zeta_prop    = {};          %vol of the vol (variance of xt)
    end
  
    methods
        function this = PMCMC_Heston(steps_mcmc, particles)
            this            = this@ParticleMarkovChainMonteCarlo(steps_mcmc, particles);
            this.mu_prop    = zeros(1, this.steps_mcmc);
            this.theta_prop = zeros(1, this.steps_mcmc);
            this.kappa_prop = zeros(1, this.steps_mcmc);
            this.zeta_prop  = zeros(1, this.steps_mcmc);
            
        end
 
        function [log_val] = call_likelihood(this, hmm, fieldname, fieldval, step)
            switch fieldname
                case 'mu_prop'
                    log_val = BootstrapParticleFilter_Heston(hmm.y, fieldval, this.theta_prop(step-1), this.kappa_prop(step-1), this.zeta_prop(step-1), this.particles);
                case 'theta_prop'
                    log_val = BootstrapParticleFilter_Heston(hmm.y, this.mu_prop(step-1), fieldval, this.kappa_prop(step-1), this.zeta_prop(step-1), this.particles);
                case 'kappa_prop'
                    log_val = BootstrapParticleFilter_Heston(hmm.y, this.mu_prop(step-1), this.theta_prop(step-1), fieldval, this.zeta_prop(step-1), this.particles);
                case 'zeta_prop'
                    log_val = BootstrapParticleFilter_Heston(hmm.y, this.mu_prop(step-1), this.theta_prop(step-1), this.kappa_prop(step-1), fieldval, this.particles);
            end
        end

        function [log_val] = call_likelihood_denom(this, hmm, step)
           log_val = BootstrapParticleFilter_Heston(hmm.y, this.mu_prop(step-1), this.theta_prop(step-1), this.kappa_prop(step-1), this.zeta_prop(step-1), this.particles); 
        end

        
        function [val] = call_prior(~, fieldname)
            switch fieldname
                case 'mu_prop'
                    val = randn*0.05;
                    %val = rand*10;
                case 'theta_prop'
                    val = rand;
                    %val = 1/gamrnd(1, 1/1); %inverse gamma
                case 'kappa_prop'
                    val = rand;
                case 'zeta_prop'
                    val = rand;
                    %val = 1/gamrnd(1, 1/1); %inverse gamma
            end
        end
        
    end
end