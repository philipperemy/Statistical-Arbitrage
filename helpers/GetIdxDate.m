function [ techid ] = GetIdxDate( pp, date_str )
    
    for dateid = pp.dt'
        if(strcmp(datestr(dateid), date_str))
            break;
        end
    end
    techid = find(dateid == pp.dt);
end

