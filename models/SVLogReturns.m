classdef SVLogReturns < HiddenMarkovModel
    properties(GetAccess = 'public', SetAccess = 'private')
        steps = 0;     %Number of steps
        beta  = 0;
    end
  
    methods

        function this = SVLogReturns(log_returns, beta)
            this = this@HiddenMarkovModel(zeros(1,length(log_returns)), log_returns);
            this.steps = zeros(1,length(log_returns));
            this.beta = beta;
        end
    end
end