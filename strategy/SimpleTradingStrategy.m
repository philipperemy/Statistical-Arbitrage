function [ pl, balance_cum ] = SimpleTradingStrategy( pp, Spread, start_idx, end_idx, boll_conf, spread_vol, disp, considerTrdCost)
    
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
    
    [ pl, balance_cum ] = Strategy_Simulator( pp, beg, T, 10000, Spread, spr, spr_uppr_trends, spr_lowr_trends, spr_lowr_trends, spr_uppr_trends, disp, considerTrdCost );

end