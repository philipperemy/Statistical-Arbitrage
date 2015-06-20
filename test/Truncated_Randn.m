function [ val ] = Truncated_Randn( mean, std, min, max )
    
    val = normrnd(mean, std);
    if(val < min)
        val = min;
    end
    
    if(max < val)
       val = max; 
    end
end

