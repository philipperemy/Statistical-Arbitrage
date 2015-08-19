function [ portfolio_cumsum ] = Cumsum_PL( pp, spreads, start_idx, end_idx, spread_ids, zscore_conf )
    portfolio_cumsum = zeros(1, length(spreads(spread_ids(1)).px));
    for id = spread_ids
        [~,cumsum_pl] = SimpleTradingStrategyZScore( pp, spreads(id), 1, 0, zscore_conf, 0, 1 );
         portfolio_cumsum = portfolio_cumsum + cumsum_pl;
    end
    portfolio_cumsum = portfolio_cumsum /length(spread_ids);
    
    if(end_idx == 0)
       end_idx = length(portfolio_cumsum);
    end
    
    portfolio_cumsum = portfolio_cumsum(start_idx:end_idx);
    portfolio_cumsum = cumsum(diff(portfolio_cumsum)) + 10000;
end

