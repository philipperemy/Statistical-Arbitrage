function [ pl, balance_cum ] = SimpleTradingStrategy( pp, Spread, start_idx, end_idx, boll_conf, disp)
    
    try
        spr = Spread.px;
    catch
        spr = Spread;
    end
    
    spr = spr(start_idx:end_idx);
    T = length(spr);
    
	[mid, uppr, lowr] = bollinger(spr, boll_conf.wsize, boll_conf.wts, boll_conf.nstd);
    
    beg = 100;
    [ spr_uppr_trends ] = CrossingPointsCalculator(spr - uppr, beg, T);
    [ spr_mid_trends  ] = CrossingPointsCalculator(spr - mid , beg, T);
    [ spr_lowr_trends ] = CrossingPointsCalculator(spr - lowr, beg, T);
    
    UP 					= 1;
    DOWN 				= -1;
	
	BUY					= 1;
	SELL				= -1;
	
    pos_buy 			= false;
    pos_sell 			= false;
    
    last_px_sell 		= 0;
    last_px_buy 		= 0;
    
    pl 					= zeros(1,T);
    
    balance_init 		= 10000;
    balance 			= balance_init;
	
	max_ss_per			= 0.5; %max short selling percentage of outstanding capital
    
    %BUG: if two positions are open at the same time, the balance is like
    %20000, more than the 10000 initial.
	for i = beg:T
        
        balance = balance + pl(i-1);
        
        % Strategy is:
        % SELL-1  : DOWN between Uppr and Spr
        % SELL-2  : DOWN between Spr  and Low  (unwind)
        % BUY-1   : UP between   Lowr and Spr
        % BUY-2   : UP between   Spr and  Uppr (unwind)
        % Only two positions BUY/SELL can be initiated at a time
        
        if(spr_uppr_trends(i) == DOWN && ~pos_sell)
            if(disp) fprintf('[%i] Initiating SELL at price %f\n', i, spr(i)); end;
            order_sell = SpreadBuildOrder( pp, Spread, SELL, balance, max_ss_per, i, disp );
            pos_sell = true; %init pos
            last_px_sell = spr(i);
			continue;
        end
        
        if(spr_lowr_trends(i) == DOWN && pos_sell)
			pos_sell = false; %unwind pos
			pl(i) = (last_px_sell - spr(i)) * order_sell.spr_qty;
			if(disp) fprintf('[%i] unwind SELL at price %f, diff %f\n', i, spr(i), pl(i)); end;
			continue;
        end
        
        if(spr_lowr_trends(i) == UP && ~pos_buy)
			if(disp) fprintf('[%i] Initiating BUY at price %f\n', i, spr(i)); end;
			order_buy = SpreadBuildOrder( pp, Spread, BUY, balance, max_ss_per, i, disp );
			pos_buy = true;
			last_px_buy = spr(i);
			continue;
        end

        if(spr_uppr_trends(i) == UP && pos_buy)
			pos_buy = false;
			pl(i) = (spr(i) - last_px_buy) * order_buy.spr_qty;
			if(disp) fprintf('[%i] unwind BUY at price %f diff %f\n', i, spr(i), pl(i)); end;
			continue;
        end
    end
    
    balance_cum = cumsum(pl) + balance_init;

end