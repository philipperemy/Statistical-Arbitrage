classdef HestonModel < HiddenMarkovModel
    properties(GetAccess = 'public', SetAccess = 'private')
        mu      = 0.91;     %rate of return of the asset 
        theta   = 1;        %long variance
        kappa   = 0.1;      %rate at which xt reverts to theta
        zeta    = 1;        %vol of the vol (variance of xt)
        steps   = 1000;     %number of steps
    end
  
    methods

        function this = HestonModel(mu, theta, kappa, zeta, steps)
            this = this@HiddenMarkovModel();
            if nargin == 5
                this.mu 	= mu;
                this.theta  = theta;
                this.kappa  = kappa;
                this.zeta   = zeta;
                this.steps  = steps;
            end
            generate(this);
        end
 
        function generate(this)
            x = zeros(1, this.steps+1);
            y = zeros(1, this.steps);
            x(1) = this.theta;      

            for i=1:this.steps
                x(i+1) = x(i) + this.kappa*(this.theta - x(i)) + this.zeta*sqrt(x(i))*randn;
                y(i)   = this.mu + sqrt(x(i))*randn;
                fprintf('x=%f,y=%f\n', x(i), y(i));
            end
            this.x = x(1:this.steps); 
            this.y = y;
        end
    end
end