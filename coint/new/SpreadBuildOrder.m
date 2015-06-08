function [ order ] = SpreadBuildOrder( pp, spread, direction, nominal, max_ss_per, t )
    %SpreadSendOrder(pp, [713, 949, 1187], [0.628830509030226, 0.148724199923458], 1, 1000, 0.5, 10)	
	
    BUY = 1;
    access_price = @(x, varargin) x(varargin{:});
    
    short_selling_nominal_limit = max_ss_per * nominal;
    
    s1_px = access_price(GetPrice(pp, spread.tuple(1)), t);
    s2_px = access_price(GetPrice(pp, spread.tuple(2)), t);
    s3_px = access_price(GetPrice(pp, spread.tuple(3)), t);
    
    beta = spread.beta;
    beta(2) = -beta(2); %easier to proceed this way
    beta(3) = -beta(3);
    
    if(direction == BUY)
        buying_nominal = s1_px*1;
        borrowing_nominal = s2_px*beta(2) + s3_px*beta(3);
        spr_qty = 1;
        while(buying_nominal < nominal && borrowing_nominal < short_selling_nominal_limit)
            spr_qty = spr_qty + 1;
            buying_nominal = s1_px*1 * spr_qty;
            borrowing_nominal = (s2_px*beta(2) + s3_px*beta(3)) * spr_qty;
        end
        spr_qty = spr_qty - 1;
        fprintf('BUY %i Spr, buying_nominal = %f, borrowing_nominal = %f, balance = %f\n', spr_qty, buying_nominal, borrowing_nominal, nominal);
        
    else
        borrowing_nominal = s1_px*1;
        buying_nominal = s2_px*beta(2) + s3_px*beta(3);
        spr_qty = 1;
        while(buying_nominal < nominal && borrowing_nominal < short_selling_nominal_limit)
            spr_qty = spr_qty + 1;
            borrowing_nominal = s1_px*1 * spr_qty;
            buying_nominal = (s2_px*beta(2) + s3_px*beta(3)) * spr_qty;
        end
        spr_qty = spr_qty - 1;
        fprintf('SELL %i Spr, borrowing_nominal = %f, buying_nominal = %f, balance = %f\n', spr_qty, borrowing_nominal, buying_nominal, nominal);
    end
	
	order = struct('spr_qty', spr_qty, 'direction', direction, 'beta', beta);

end


% 1,-2,4
% 
% 100; 140; 60
% 
% BUY SPR@60
% 1*100 - 2*140 + 60*4 = 60
% BUY 1-1@100
% SELL 2-2@140
% BUY 3-4@60
% 
% Cash involved: 340 EUR
% Cash borrowed: 280 EUR
% 
% Q*(1*100-2*140+60*4) = Q*(1*100-beta(2)*140+beta(3)*60)
% Q = quantity to buy or sell
% Q*(340-280)
% 
% 100
