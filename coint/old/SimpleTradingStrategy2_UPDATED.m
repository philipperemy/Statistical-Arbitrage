function [ pl ] = SimpleTradingStrategy2_UPDATED( Spread, start_idx, end_idx, boll_conf)
    
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    spr = spr(start_idx:end_idx);
    [mid, uppr, lowr] = bollinger(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd);
    plot([mid, uppr, lowr, spr]);
    
    T = length(spr);
    beg = 150;
    [ spr_uppr_trends ] = CrossingTrends(spr, uppr, beg, T);
    [ spr_lowr_trends ] = CrossingTrends(spr, lowr, beg, T);
    
    UP = 1;
    DOWN = -1;
    pos_buy = false;
    pos_sell = false;
    last_px = 0;
    pl = zeros(1,T);

    for i = beg:T
        % Strategy is:
        % SELL-1  : DOWN between Uppr and Spr
        % SELL-2  : DOWN between Low and Spr
        % BUY-1   : UP between Lowr and Spr
        % BUY-2   : UP between Uppr and Spr
        % When a position is initiated, no other positions are allowed
        
        if(spr_uppr_trends(i) == DOWN && ~pos_buy && ~pos_sell)
            pos_sell = true; %init pos
            last_px = spr(i);
            continue;
        end
        
        if(spr_lowr_trends(i) == DOWN && pos_sell && ~pos_buy)
           pos_sell = false; %unwind pos
           pl(i) = last_px - spr(i);
           continue;
        end
        
        if(spr_lowr_trends(i) == UP && ~pos_buy && ~pos_sell)
           pos_buy = true;
           last_px = spr(i);
           continue;
        end
        
        if(spr_uppr_trends(i) == UP && pos_buy && ~pos_sell)
           pos_buy = false;
           pl(i) = spr(i) - last_px;
           continue;
        end
    end

end