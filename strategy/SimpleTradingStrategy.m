function [ pl, balance_cum, trds ] = SimpleTradingStrategy( pp, Spread, start_idx, end_idx, boll_conf, spread_vol, disp, considerTrdCost, strat_simulator, initial_bet)
    
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    if(end_idx ~= 0)
        spr = spr(start_idx:end_idx);
    end
    
    T = length(spr);
    
    if(spread_vol ~= 0)
        [~, uppr, lowr] = Compute_Bollinger_Bands_SV(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd, spread_vol);
    else
        [~, uppr, lowr] = bollinger(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd);
    end
	
    beg = boll_conf.wsize;
    [ spr_uppr_trends ] = CrossingPointsCalculator(spr - uppr, beg, T);
    %[ spr_mid_trends  ] = CrossingPointsCalculator(spr - mid , beg, T);
    [ spr_lowr_trends ] = CrossingPointsCalculator(spr - lowr, beg, T);
    
    [ pl, balance_cum, trds ] = strat_simulator( pp, beg, T, initial_bet, Spread, spr, spr_uppr_trends, spr_lowr_trends, spr_lowr_trends, spr_uppr_trends, disp, considerTrdCost );

    %stopping point
    if(any(balance_cum < 0))
        id = find( balance_cum < 0); id = id(1);
        balance_cum(id:end) = zeros(1,length(balance_cum)-id+1);
    end
    
end