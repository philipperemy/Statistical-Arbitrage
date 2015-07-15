function [ pp_train, pp_test ] = Create_Sets( pp, dateIdx_start, dateIdx_end, percentage )
    if(percentage > 1)
        percentage = percentage/100;
    end

    dateIdx_split = round(dateIdx_start + (dateIdx_end-dateIdx_start)*percentage);
    
    pp_train = TruncateData(pp, dateIdx_start, dateIdx_split);
    pp_test  = TruncateData(pp, dateIdx_split+1, dateIdx_end);
end

