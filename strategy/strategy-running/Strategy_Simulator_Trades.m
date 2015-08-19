function [ pl, balance_cum ] = Strategy_Simulator_Trades( pp, beg, T, balance_init, Spread, spr, sell_open, sell_close, buy_open, buy_close, disp, considerTrdCost )

    UP 					= 1;
    DOWN 				= -1;
	
	BUY					= 1;
	SELL				= -1;
	
    pos_buy 			= false;
    pos_sell 			= false;
    
    last_px_sell 		= 0;
    last_px_buy 		= 0;
    
    pl 					= zeros(1,T);
    balance 			= balance_init;
	max_ss_per			= 0.5; %max short selling percentage of outstanding capital
    
    for i = beg:T
        
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
			pl(i) = ( (last_px_sell - spr(i)) * order_sell.spr_qty ) - transactionCost;
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
			pos_buy = false;
            
            if(considerTrdCost == 1)
                transactionCost = ( order_buy.spr_qty * last_px_buy * (5/10000) );
            else
                transactionCost = 0;
            end
            
			pl(i) = ( (spr(i) - last_px_buy) * order_buy.spr_qty ) - transactionCost;
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

