function [ mat ] = GetPriceArray( pp, ids )
    mat = [];
    
    if(size(ids, 2) == 1)
       ids = ids'; 
    end
    
    for i = ids
        mat = [mat, GetPrice(pp, i)];
    end
end

