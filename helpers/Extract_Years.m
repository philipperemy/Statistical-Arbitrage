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
        seq_twobytwo = [1 2*(1:1:floor(length(techids)/2))+1];
        if(seq_twobytwo(end) > length(techids))
            seq_twobytwo(end) = [];
        end
       techids = techids(seq_twobytwo); %select period of two years 
    end
end

