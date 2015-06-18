classdef StochasticVolatilityModel < HiddenMarkovModel
    properties(GetAccess = 'public', SetAccess = 'private')
        rho = 0.91;       %Scale parameter for X
        sigma = 1;        %Standard deviation of the X (process)
        beta = 0.5;       %Scale parameter for Y
        steps = 1000;     %Number of steps
    end
  
    methods

        function this = StochasticVolatilityModel(rho, sigma, beta, steps, nu)
            this = this@HiddenMarkovModel();
            if nargin == 5
                this.rho = rho;
                this.sigma = sigma;
                this.beta = beta;
                this.steps = steps;
                this.nu = nu;
            end
            generate(this);
        end
 
        function generate(this)
            x = zeros(1, this.steps+1);
            y = zeros(1, this.steps);
            x(1) = randn * sqrt(this.sigma^2/(1-this.rho^2));      

            for i=1:this.steps
                x(i+1) = this.rho  * x(i)  + this.sigma*randn;
                y(i)   = this.beta * exp(0.5*x(i))*randn;
            end
            this.x = x(1:this.steps); 
            this.y = y;
        end
    end
end