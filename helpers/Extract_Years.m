%Extract the technical ids of the first day of each year
function [ techids ] = Extract_Years(pp, twobytwo)
    %Be careful of Leap year
    techids = [];
    for dateid = pp.dt'
        if(strstartswith(datestr(dateid), '01-Jan'))
            techids(end+1) = find(dateid == pp.dt);
        end
    end
    
    if(exist('twobytwo', 'var') && twobytwo)
       techids = techids([1 2*(1:1:floor(length(techids)/2))+1]); %select period of two years 
    end
end

