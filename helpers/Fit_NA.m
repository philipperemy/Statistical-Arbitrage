function [ vec ] = Fit_NA( vec )
    
    N = length(vec);
    for i = 2:(N-1)
       if(isnan(vec(i)) || vec(i) == 0)
          vec(i) = 0.5 * vec(i-1) + 0.5 * vec(i+1);
       end
    end

end

