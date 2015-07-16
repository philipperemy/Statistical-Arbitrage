%Extract the technical ids of the first day of each year
function [ techids ] = Extract_Years(pp)
    %Be careful of Leap year
    techids = [];
    for dateid = pp.dt'
        if(strstartswith(datestr(dateid), '01-Jan'))
            techids(end+1) = find(dateid == pp.dt);
        end
    end
end

