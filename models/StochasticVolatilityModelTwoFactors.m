classdef StochasticVolatilityModelTwoFactors < HiddenMarkovModel
    properties(GetAccess = 'public', SetAccess = 'private')
        rho = 0.91;       %Scale parameter for X
        sigma = 1;        %Standard deviation of the X (process)
        
        rho2 = 0.5;
        sigma2 = 3;
        
        beta = 0.5;       %Scale parameter for Y
        steps = 2000;     %Number of steps
        z;
    end
  
    methods

        function this = StochasticVolatilityModelTwoFactors()
            this = this@HiddenMarkovModel();
            generate(this);
        end
 
        function generate(this)
            x = zeros(1, this.steps+1);
            x(1) = randn * sqrt(this.sigma^2/(1-this.rho^2));      
            
            this.z = zeros(1, this.steps+1);
            this.z(1) = randn * sqrt(this.sigma2^2/(1-this.rho2^2));  
            
            y = zeros(1, this.steps);
            for i=1:this.steps
                x(i+1)      = this.rho  * x(i)       + this.sigma*randn;
                this.z(i+1) = this.rho2 * this.z(i)  + this.sigma2*randn;
                
                y(i)        = this.beta * exp(0.5*(x(i)+this.z(i)))* randn;
            end
            this.x = x(1:this.steps); 
            this.z = this.z(1:this.steps); 
            this.y = y;
        end
    end
end