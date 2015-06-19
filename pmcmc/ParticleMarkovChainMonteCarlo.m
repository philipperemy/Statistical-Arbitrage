classdef ParticleMarkovChainMonteCarlo < handle
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        steps_mcmc = 100;
        particles = 1000;
        
        NO_UPDATE = 0;
        UPDATE = 1;
    end
	
	methods (Abstract)
		log_val = call_likelihood(this, hmm, fieldname, fieldval, step);
        log_val = call_likelihood_denom(this, hmm, step);
		val = call_prior(this, fieldname);
	end
  
    methods
        function this = ParticleMarkovChainMonteCarlo(steps_mcmc, particles)
            if nargin == 2
                this.steps_mcmc = steps_mcmc; 
                this.particles = particles;
            end
        end
 
        function [ret] = check_update(this, log_numerator, log_denominator)
            ret = this.NO_UPDATE;
            if (rand < min(1, exp(log_numerator - log_denominator)))
                ret = this.UPDATE;
            end
        end
        
        function this = update_field(this, log_numerator, log_denominator, step, fieldname, value)
            new_val = value;
            if(check_update(this, log_numerator, log_denominator) == this.NO_UPDATE)
                new_val = this.(fieldname)(step - 1);
            end
            
            this.(fieldname)(step) = new_val;
        end

        %Can be used as the base in the class
        %make call_likelihood abstract
        function this = run(this, hmm)
           
            i = 1;
            fields = {};
            for field = fieldnames(this)'
                fieldname = char(field);
                if(strendswith(fieldname, '_prop'))
                    fields(i) = cellstr(fieldname);
                    i = i + 1;
                    
                    %init field with prior
                    this.(fieldname)(1) = this.call_prior(fieldname);
                end
            end
            
            for step = 2:this.steps_mcmc

                %Core
                log_denominator = this.call_likelihood_denom(hmm, step);
                for i = 1:length(fields)
                    fieldname = fields{i};
                    if(strendswith(fieldname, '_prop'))
                        
                        while true
                            prior_val = this.call_prior(fieldname);
                            log_numerator = this.call_likelihood(hmm, fieldname, prior_val, step);
                            
                            if(~isnan(log_numerator) && log_numerator ~= -Inf)
                               break; 
                            end
                          
                        end
                        
                        this.update_field(log_numerator, log_denominator, step, fieldname, prior_val);
                    end
                end

                %print
                if(mod(step, 1) == 0)
                    str = 'MCMC step: %i , ';
                    vals = [step];
                    for i = 1:length(fields)
                        fieldname = fields{i};
                        ref = this.(fieldname);
                        val = ref(step-1);
                        vals = [vals, val];
                        str = strcat(str, ' ', fieldname, ' = %f , ');
                    end
                    str = strcat(str, ' likelihood = %f');
                    str = strcat(str, '\n');
                    vals = [vals, log_denominator];
                    fprintf(str, vals);
                end
            end
        end
    end
end