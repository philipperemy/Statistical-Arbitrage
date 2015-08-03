function [ stocks, stocks_tst ] = Prepare_Set( stocks, stocks_tst )

    %synchronise training and testing sets

    [~,unique_idx]      =   unique(stocks.px_clean','rows');
    stocks              =   TruncateDataIds( stocks, unique_idx );
    stocks_tst          =   TruncateDataIds( stocks_tst, unique_idx );

    stocks_count = size(stocks.px,2);
    lens = zeros(1,stocks_count);
    for i = 1:stocks_count
        lens(i) = length(GetPrice(stocks, i));
    end
    [no,xo]  = hist(lens,unique(lens));
    [~,I]   = sort(no, 'descend');
    vals    = xo(I);
    if(vals(1) == 0)
       mode_val = vals(2);
    else
       mode_val = vals(1);
    end

    stocks_idx          =   find(lens == mode_val);
    
    stocks              =   TruncateDataIds( stocks, stocks_idx );
    stocks_tst          =   TruncateDataIds( stocks_tst, stocks_idx );

    %precached returns - stationary time series
    stocks.ret          = zeros(mode_val-1, length(stocks_idx));
    for i = 1:length(stocks_idx)
        stocks.ret(:,i) = diff(log(GetPrice(stocks, i)));
    end
end

