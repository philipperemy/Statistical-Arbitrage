function [ pp ] = TruncateDataIds( pp, ids )

    pp.sector       =   pp.sector(ids);    
    pp.tickers      =   pp.tickers(ids);
    if(~isempty(pp.w))
        pp.w        =   pp.w(:,ids);
    end
    pp.names        =   pp.names(ids);
    pp.px           =   pp.px(:,ids);

end