function [ portfolio_cumsum ] = Cumsum_PL( pp, spreads, spread_ids, zscore_conf )
    portfolio_cumsum = zeros(1, length(spreads(spread_ids(1)).px));
    for id = spread_ids
        [~,cumsum_pl] = SimpleTradingStrategyZScore( pp, spreads(id), 1, 0, zscore_conf, 0, 1 );
         portfolio_cumsum = portfolio_cumsum + cumsum_pl;
    end
    portfolio_cumsum = portfolio_cumsum /length(spread_ids);

end

