classdef HiddenMarkovModel < handle
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        x; %process
        y; %measurement
    end
  
    methods

        function this = HiddenMarkovModel(x, y)
            if nargin == 2
                this.x = x;
                this.y = y;
            end
        end
 
    end
end