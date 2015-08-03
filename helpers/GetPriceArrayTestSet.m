function [ mat_ret ] = GetPriceArrayTestSet( pp, ids )
    mat = [];

    if(size(ids, 2) == 1)
       ids = ids'; 
    end

    for i = ids
        mat = [mat, pp.px(:,i)];
    end
    
    mat_ret = [];
    for i = 1:length(mat)
        if(~any(isnan(mat(i,:))))
            mat_ret = [mat_ret; mat(i,:)];
        end
    end
end

