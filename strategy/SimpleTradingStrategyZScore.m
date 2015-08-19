function [ pl, balance_cum ] = SimpleTradingStrategyZScore( pp, Spread, start_idx, end_idx, zscore_conf, disp, considerTrdCost)
    
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
    
    [ sell_open ] = CrossingPointsCalculator(spr_n - zscore_conf.sell_open, 2, T);
    [ sell_close  ] = CrossingPointsCalculator(spr_n - zscore_conf.sell_close , 2, T);
    [ buy_open ] = CrossingPointsCalculator(spr_n - zscore_conf.buy_open, 2, T);
    [ buy_close  ] = CrossingPointsCalculator(spr_n - zscore_conf.buy_close , 2, T);
    
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
    for i = 2:T
        
        balance = balance + pl(i-1);
        
        if(sell_open(i) == DOWN && ~pos_sell)
            if(disp) fprintf('[%i] Initiating SELL at price %f\n', i, spr(i)); end;
            order_sell = SpreadBuildOrder( pp, Spread, SELL, balance, max_ss_per, i, disp );
            pos_sell = true; %init pos
            last_px_sell = spr(i);
			continue;
        end
        
        if(sell_close(i) == DOWN && pos_sell)
            
            if(considerTrdCost == 1)
                transactionCost = ( order_sell.spr_qty * last_px_sell * (5/10000) );
            else
                transactionCost = 0;
            end
            
			pos_sell = false; %unwind pos
			pl(i) = ((last_px_sell - spr(i)) * order_sell.spr_qty) - transactionCost;
			if(disp) fprintf('[%i] unwind SELL at price %f, diff %f\n', i, spr(i), pl(i)); end;
			continue;
        end
        
        if(buy_open(i) == UP && ~pos_buy)
			if(disp) fprintf('[%i] Initiating BUY at price %f\n', i, spr(i)); end;
			order_buy = SpreadBuildOrder( pp, Spread, BUY, balance, max_ss_per, i, disp );
			pos_buy = true;
			last_px_buy = spr(i);
			continue;
        end

        if(buy_close(i) == UP && pos_buy)
            
            if(considerTrdCost == 1)
                transactionCost = ( order_buy.spr_qty * last_px_buy * (5/10000) );
            else
                transactionCost = 0;
            end
            
			pos_buy = false;
			pl(i) = ((spr(i) - last_px_buy) * order_buy.spr_qty) - transactionCost;
			if(disp) fprintf('[%i] unwind BUY at price %f diff %f\n', i, spr(i), pl(i)); end;
			continue;
        end
    end
    
    %close position at the end
    if(pos_buy)
        if(considerTrdCost == 1)
            transactionCost = ( order_buy.spr_qty * last_px_buy * (5/10000) );
        else
            transactionCost = 0;
        end
        
        pl(i) = ((spr(i) - last_px_buy) * order_buy.spr_qty) - transactionCost;
        if(disp) fprintf('[%i] unwind BUY at price %f diff %f\n', i, spr(i), pl(i)); end;
    end
    
    if(pos_sell)
        if(considerTrdCost == 1)
            transactionCost = ( order_sell.spr_qty * last_px_sell * (5/10000) );
        else
            transactionCost = 0;
        end

        pl(i) = ((last_px_sell - spr(i)) * order_sell.spr_qty) - transactionCost;
        if(disp) fprintf('[%i] unwind SELL at price %f, diff %f\n', i, spr(i), pl(i)); end;
    end
    
    balance_cum = cumsum(pl) + balance_init;
end