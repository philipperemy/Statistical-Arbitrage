%Keep only rows between start and end
function [ pp_ret ] = RetrieveStocksSector( pp, sector_name )
    ids_sector          =   strmatch(sector_name, pp.sector);
    pp_ret              =   pp;
    pp_ret.sector       =   pp.sector(ids_sector);    
    pp_ret.tickers      =   pp.tickers(ids_sector);
    if(isfield('pp', 'w') && ~isempty(pp.w))
        pp_ret.w        =   pp.w(:,ids_sector);
    end
    pp_ret.names        =   pp.names(ids_sector);
    pp_ret.px           =   pp.px(:,ids_sector);
    pp_ret.ret          =   pp.ret(:,ids_sector);
end

