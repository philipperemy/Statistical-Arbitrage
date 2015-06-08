function [ spr_qty ] = SpreadSendOrder( pp, tuples, beta, direction, nominal, max_short_selling_per_nominal, t )
    %SpreadSendOrder(pp, [713, 949, 1187], [0.628830509030226, 0.148724199923458], 1, 1000, 0.5, 10)

    %Only for beta(1,2)>0
    
    short_selling_nominal_limit = max_short_selling_per_nominal*nominal;
    
    s1_t = Convert(pp, tuples(1));
    s1_t = s1_t(t);
    
    s2_t = Convert(pp, tuples(2));
    s2_t = s2_t(t);
    
    s3_t = Convert(pp, tuples(3));
    s3_t = s3_t(t);
    
    if(direction == 1)
        buying_nominal = s1_t*1;
        borrowing_nominal = s2_t*beta(1) + s3_t*beta(2);
        spr_qty = 1;
        while(buying_nominal < nominal && borrowing_nominal < short_selling_nominal_limit)
            spr_qty = spr_qty + 1;
            buying_nominal = s1_t*1 * spr_qty;
            borrowing_nominal = (s2_t*beta(1) + s3_t*beta(2)) * spr_qty;
        end
        spr_qty = spr_qty - 1;
        fprintf('BUY %i Spr, buying_nominal = %f, borrowing_nominal = %f, balance = %f\n', spr_qty, buying_nominal, borrowing_nominal, nominal);
        
    else
        borrowing_nominal = s1_t*1;
        buying_nominal = s2_t*beta(1) + s3_t*beta(2);
        spr_qty = 1;
        while(buying_nominal < nominal && borrowing_nominal < short_selling_nominal_limit)
            spr_qty = spr_qty + 1;
            borrowing_nominal = s1_t*1 * spr_qty;
            buying_nominal = (s2_t*beta(1) + s3_t*beta(2)) * spr_qty;
        end
        spr_qty = spr_qty - 1;
        fprintf('SELL %i Spr, borrowing_nominal = %f, buying_nominal = %f, balance = %f\n', spr_qty, borrowing_nominal, buying_nominal, nominal);
    end

end




% 1,-2,4
% 
% 1000
% 
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
% Q*(1*100-2*140+60*4) = Q*(1*100-beta(1)*140+beta(2)*60)
% Q = quantity to buy or sell
% Q*(340-280)
% 
% 100
