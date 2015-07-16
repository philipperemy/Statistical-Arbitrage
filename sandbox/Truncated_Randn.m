function [ val ] = Truncated_Randn( mean, std, min, max )
    
    val = normrnd(mean, std);
    
    if(min == 0)
       min = 1e-5; 
    end
    
    if(val < min)
        val = min;
    end
    
    if(max < val)
       val = max; 
    end
end

