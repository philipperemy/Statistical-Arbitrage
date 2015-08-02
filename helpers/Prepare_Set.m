function [ stocks ] = Prepare_Set( stocks )

    [~,unique_idx]      =   unique(stocks.px_clean','rows');
    stocks.sector       =   stocks.sector(unique_idx);    
    stocks.tickers      =   stocks.tickers(unique_idx);
    
    if(~isempty(stocks.w))
        stocks.w        =   stocks.w(:,unique_idx);
    end

    stocks.names        =   stocks.names(unique_idx);
    stocks.px           =   stocks.px(:,unique_idx);

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
    stocks.sector       =   stocks.sector(stocks_idx);
    stocks.tickers      =   stocks.tickers(stocks_idx);
    stocks.names        =   stocks.names(stocks_idx);
    stocks.px           =   stocks.px(:,stocks_idx);
    if(~isempty(stocks.w))
        stocks.w        =   stocks.w(:,stocks_idx);
    end
    
    %precached returns - stationary time series
    stocks.ret          = zeros(mode_val-1, length(stocks_idx));
    for i = 1:length(stocks_idx)
        stocks.ret(:,i) = diff(log(GetPrice(stocks, i)));
    end
end

