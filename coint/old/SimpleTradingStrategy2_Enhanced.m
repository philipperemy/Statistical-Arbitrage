%http://quant.stackexchange.com/questions/17295/limits-on-short-selling
%Don't shortsell more than %X% percent of your total outstanding nominal

function [ pl, balance ] = SimpleTradingStrategy2_Enhanced( pp, tuples, start_idx, end_idx, boll_conf, initial_cash)
    
    [ Spread, beta, tuples_ret ] = SpreadConstructor(pp, tuples);
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    spr = spr(start_idx:end_idx);
    [mid, uppr, lowr] = bollinger(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd);
    
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

    max_short_selling_per_nominal = 0.5;
    qty = 0;
    SELL = 0;
    BUY = 1;
    
    current_cash = initial_cash;
    
    for i = beg:T
        
        balance(i) = current_cash;
        current_cash = current_cash + pl(i-1);
        
        % Strategy is:
        % SELL-1  : DOWN between Uppr and Spr
        % SELL-2  : DOWN between Low and Spr
        % BUY-1   : UP between Lowr and Spr
        % BUY-2   : UP between Uppr and Spr
        % When a position is initiated, no other positions are allowed
        
        if(spr_uppr_trends(i) == DOWN && ~pos_buy && ~pos_sell)
            pos_sell = true; %init pos - OPEN SELL
            qty = SpreadSendOrder( pp, tuples_ret, beta, SELL, current_cash, max_short_selling_per_nominal, i );
            last_px = spr(i);
            continue;
        end
        
        if(spr_lowr_trends(i) == DOWN && pos_sell && ~pos_buy)
           pos_sell = false; %unwind pos - CLOSE SELL
           pl(i) = (last_px - spr(i))*qty;
           qty = 0;
           continue;
        end
        
        if(spr_lowr_trends(i) == UP && ~pos_buy && ~pos_sell)
           pos_buy = true; %init pos - OPEN BUY
           qty = SpreadSendOrder( pp, tuples_ret, beta, BUY, current_cash, max_short_selling_per_nominal, i );
           last_px = spr(i); %last value is not exactly right. Careful
           continue;
        end
        
        if(spr_uppr_trends(i) == UP && pos_buy && ~pos_sell)
           pos_buy = false; %unwind pos - CLOSE BUY
           pl(i) = (spr(i) - last_px)*qty;
           qty = 0;
           continue;
        end
    end

end