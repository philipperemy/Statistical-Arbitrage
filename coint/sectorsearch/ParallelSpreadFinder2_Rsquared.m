run('../Init.m');
addpath('../../helpers');
diary off;
delete 'log.txt';
diary 'log.txt';
tic;

year_start_dates    = Extract_Years(pp);
year_2013_start     = year_start_dates(end-1);
year_2013_end       = year_start_dates(end);

[pp_train_2013, pp_test_2013] = Create_Sets( pp, year_2013_start, year_2013_end, 0.7 );

%do not use px_clean! it's biased.
stocks          =   pp_train_2013;
spreads_count   =   1;
sectors         =   unique(pp.sector);

%
% For now R squared computed on stocks directly and not on stationary
% series
stocks.px = diff(stocks.px);

%count total more than 2billion
% s_count = 1;
% length_ids = 108;
% for i = 1:length_ids
%     for j = 1:length_ids
%         for k = (j+1):length_ids
%             for l = (k+1):length_ids
%                 for m = (l+1):length_ids
%                     s_count = s_count + 1;
%                 end
%             end
%         end
%     end
% end

%Financials
for sector = sectors(4)' %REMOVE THIS SHIT
    ids = strmatch(sector, pp.sector);
    
    %useful to avoid going from 1:end but start from the middle
    ids_length = length(ids);
    for i = 1:ids_length
        fprintf('i = %i\n', i);
        stock_1 = GetPrice(stocks, ids(i));
        
        %[25 26 27]
        %[26 25 27] should be here
        %[26 27 25] we don't care
        %[26 26 25] impossible. move to the next one
        %So i should go through the range and j too with i diff from j
        for j = 1:ids_length 
            
            if(i == j)
                continue;
            end
            
            stock_2 = GetPrice(stocks, ids(j));
            
            for k = (j+1):ids_length
                stock_3 = GetPrice(stocks, ids(k));
                
                for l = (k+1):ids_length
                    stock_4 = GetPrice(stocks, ids(l));
                    
                    for m = (l+1):ids_length
                        stock_5 = GetPrice(stocks, ids(m));
                        
                        cor_mat = corr([stock_1 stock_2 stock_3 stock_4 stock_5]);
                        
                        if(any(size(cor_mat) ~= [5 5]))
                           continue; 
                           %that means one vector contains only NaN
                        end
                        
                        c = cor_mat(2:5,1);
                        R = corr([stock_2 stock_3 stock_4 stock_5]);
                        val = c'*inv(R)*c;
                        spreads_count = spreads_count + 1;
                        fprintf('%i, %i, %i, %i, %i, %i, %f\n', spreads_count, ids(i), ids(j), ids(k), ids(l), ids(m), val);
                    end
                end
            end
        end
    end
end

diary off;
diary on;

% 
% [~,unique_idx]  =   unique(pp.px_clean','rows');
% pp.px_clean     =   pp.px_clean(:,unique_idx);
% pp.sector       =   pp.sector(unique_idx);    
% pp.tickers      =   pp.tickers(unique_idx);
% pp.names        =   pp.names(unique_idx);
% pp.px           =   pp.px(:,unique_idx);
% pp.w            =   pp.w(:,unique_idx);
% 
% pp.px_clean(1,:)=   [];
% 
% pp.px(1,:)      =   [];
% pp.dt(1)        =   [];

