function [ pl, balance_cum ] = SimpleTradingStrategy( pp, Spread, start_idx, end_idx, boll_conf)
    
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    spr = spr(start_idx:end_idx);
    [mid, uppr, lowr] = bollinger(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd);
    %plot([mid, uppr, lowr, spr]);
    
    T = length(spr);
    beg = 100;
    [ spr_uppr_trends ] = CrossingPointsCalculator(spr - uppr, beg, T);
    [ spr_mid_trends ] = CrossingPointsCalculator(spr - mid, beg, T);
    [ spr_lowr_trends ] = CrossingPointsCalculator(spr - lowr, beg, T);
    
    UP = 1;
    DOWN = -1;
    pos_buy = false;
    pos_sell = false;
    
    last_px_sell = 0;
    last_px_buy = 0;
    
    pl = zeros(1,T);
    
    balance_init = 10000;
    balance = balance_init;
    for i = beg:T
        
        balance = balance + pl(i-1);
        
        % Strategy is:
        % SELL-1  : DOWN between Uppr and Spr
        % SELL-2  : DOWN between Mid and Spr
        % BUY-1   : UP between Lowr and Spr
        % BUY-2   : UP between Mid and Spr
        % When a position is initiated, no other positions are allowed
        
        if(spr_uppr_trends(i) == DOWN && ~pos_sell)
            fprintf('[%i] Initiating SELL at price %f\n', i, spr(i));
            order_sell = SpreadBuildOrder( pp, Spread, 0, balance, 0.5, i );
            pos_sell = true; %init pos
            last_px_sell = spr(i);
            continue;
        end
        
        if(spr_lowr_trends(i) == DOWN && pos_sell)
           pos_sell = false; %unwind pos
           pl(i) = (last_px_sell - spr(i)) * order_sell.spr_qty;
           fprintf('[%i] unwind SELL at price %f, diff %f\n', i, spr(i), pl(i));
           continue;
        end
        
        if(spr_lowr_trends(i) == UP && ~pos_buy)
           fprintf('[%i] Initiating BUY at price %f\n', i, spr(i));
           order_buy = SpreadBuildOrder( pp, Spread, 0, balance, 0.5, i );
           pos_buy = true;
           last_px_buy = spr(i);
           continue;
        end

        if(spr_uppr_trends(i) == UP && pos_buy)
           pos_buy = false;
           pl(i) = (spr(i) - last_px_buy) * order_buy.spr_qty;
           fprintf('[%i] unwind BUY at price %f diff %f\n', i, spr(i), pl(i));
           continue;
        end
    end
    
    balance_cum = cumsum(pl) + balance_init;

end