function [ pl, balance_cum, trds ] = Strategy_Simulator( pp, beg, T, balance_init, Spread, spr, sell_open, sell_close, buy_open, buy_close, disp, considerTrdCost )

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
    
    cash_balance        = balance;
    trds                = 0;
    for i = beg:T
        
        balance = balance + pl(i-1);
        if(pos_sell)
            pl(i) = (spr(i-1) - spr(i)) * order_sell.spr_qty;
        end
        
        if(pos_buy)
            pl(i) =  (spr(i) - spr(i-1)) * order_buy.spr_qty;
        end
        
        if(sell_open(i) == DOWN && ~pos_sell)
            if(disp) fprintf('[%i] Initiating SELL at price %f\n', i, spr(i)); end;
            order_sell = SpreadBuildOrder( pp, Spread, SELL, cash_balance, max_ss_per, i, disp );
            pos_sell = true; %init pos
            last_px_sell = spr(i);
            cash_balance = cash_balance - order_sell.spr_qty * last_px_sell;
            trds = trds + 1;
			continue;
        end
        
        if(sell_close(i) == DOWN && pos_sell)
            
            if(considerTrdCost == 1)
                transactionCost = ( order_sell.spr_qty * last_px_sell * (5/10000) );
            else
                transactionCost = 0;
            end
            
            cash_balance = cash_balance + order_sell.spr_qty * last_px_sell; %get back your bet
            cash_balance = cash_balance + (last_px_sell - spr(i))*order_sell.spr_qty - transactionCost; %add or diff
            pos_sell = false; %unwind pos
			pl(i) = pl(i) - transactionCost;
			if(disp) fprintf('[%i] unwind SELL at price %f, bought at %f\n', i, spr(i), last_px_sell); end;
			continue;
        end
        
        if(buy_open(i) == UP && ~pos_buy)
			if(disp) fprintf('[%i] Initiating BUY at price %f\n', i, spr(i)); end;
			order_buy = SpreadBuildOrder( pp, Spread, BUY, cash_balance, max_ss_per, i, disp );
			pos_buy = true;
			last_px_buy = spr(i);
            cash_balance = cash_balance - order_buy.spr_qty * last_px_buy;
            trds = trds + 1;
			continue;
        end

        if(buy_close(i) == UP && pos_buy)
			pos_buy = false;
            
            if(considerTrdCost == 1)
                transactionCost = ( order_buy.spr_qty * last_px_buy * (5/10000) );
            else
                transactionCost = 0;
            end
            
            cash_balance = cash_balance + order_buy.spr_qty * last_px_buy; %get back your bet
            cash_balance = cash_balance + (spr(i)-last_px_buy)*order_buy.spr_qty - transactionCost; %add or diff
			pl(i) = pl(i) - transactionCost;
			if(disp) fprintf('[%i] unwind BUY at price %f, bought at %f\n', i, spr(i), last_px_buy); end;
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

