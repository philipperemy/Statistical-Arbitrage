function [ pl, balance_cum ] = SimpleTradingStrategyZScore( pp, Spread, start_idx, end_idx, zscore_conf, disp, considerTrdCost, strat_simulator)
    
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    if(end_idx ~= 0)
        spr = spr(start_idx:end_idx);
    end
    
    T = length(spr);
    
    mu = mean(spr);
    sig = std(spr);
    spr_n = (spr - mu)/sig;
    
    [ sell_open ]   = CrossingPointsCalculator(spr_n - zscore_conf.sell_open, 2, T);
    [ sell_close  ] = CrossingPointsCalculator(spr_n - zscore_conf.sell_close , 2, T);
    [ buy_open ]    = CrossingPointsCalculator(spr_n - zscore_conf.buy_open, 2, T);
    [ buy_close  ]  = CrossingPointsCalculator(spr_n - zscore_conf.buy_close , 2, T);
    
    [ pl, balance_cum ] = strat_simulator( pp, 2, T, 10000, Spread, spr, sell_open, sell_close, buy_open, buy_close, disp, considerTrdCost );
end