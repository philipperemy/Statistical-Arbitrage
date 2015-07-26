classdef StochasticVolatilityModelMA1 < HiddenMarkovModel
    properties(GetAccess = 'public', SetAccess = 'private')
        rho = 0.91;       %Scale parameter for X
        sigma = 1;        %Standard deviation of the X (process)
        beta = 0.5;       %Scale parameter for Y
		psi	 = 0.8;
        steps = 1000;     %Number of steps
    end
  
    methods

        function this = StochasticVolatilityModelMA1(rho, sigma, beta, steps, psi)
            this = this@HiddenMarkovModel();
            if nargin == 5
                this.rho 	= rho;
                this.sigma 	= sigma;
                this.beta 	= beta;
                this.steps 	= steps;
				this.psi	= psi;
            end
            generate(this);
        end
 
        function generate(this)
            x = zeros(1, this.steps+1);
            y = zeros(1, this.steps);
            x(1) = randn * sqrt(this.sigma^2/(1-this.rho^2));      

			innov_yt_1 = randn;
            for i=1:this.steps
                x(i+1) 		= this.rho  * x(i)  + this.sigma*randn;
				innov_yt 	= randn;
                y(i)   		= this.beta * exp(0.5*x(i))* innov_yt + this.psi * this.beta * exp(0.5*x(i))* innov_yt_1;
				innov_yt_1 	= innov_yt;
            end
            this.x = x(1:this.steps); 
            this.y = y;
        end
    end
end