run('../Init.m');
addpath('../../helpers');
diary off;
delete 'Financials.txt';
diary 'Financials.txt';
tic;

year_start_dates    = Extract_Years(pp);
year_2013_start     = year_start_dates(end-1);
year_2013_end       = year_start_dates(end);

[pp_train_2013, pp_test_2013] = Create_Sets( pp, year_2013_start, year_2013_end, 0.7 );

%do not use px_clean! it's biased.
stocks          =   pp_train_2013;
sectors         =   unique(pp.sector);

%
% For now R squared computed on stocks directly and not on stationary
% series
stocks.px = diff(stocks.px);

[~,unique_idx]  =   unique(stocks.px_clean','rows');
stocks.px_clean     =   stocks.px_clean(:,unique_idx);
stocks.sector       =   stocks.sector(unique_idx);    
stocks.tickers      =   stocks.tickers(unique_idx);
stocks.names        =   stocks.names(unique_idx);
stocks.px           =   stocks.px(:,unique_idx);
% 
% 
% for sector = sectors(4)' %REMOVE THIS SHIT
%     ids = strmatch(sector, stocks.sector);
% 
%     i = find(ids == 65);
%     j = find(ids == 61);
%     k = find(ids == 65);
%     l = find(ids == 565);
%     m = find(ids == 566);
%     
%     stock_1 = GetPrice(stocks, ids(i));
%     stock_2 = GetPrice(stocks, ids(j));
%     stock_3 = GetPrice(stocks, ids(k));
%     stock_4 = GetPrice(stocks, ids(l));
%     stock_5 = GetPrice(stocks, ids(m));
%     cor_mat = corr([stock_1 stock_2 stock_3 stock_4 stock_5]);
% 
%     c = cor_mat(2:5,1);
%     R = corr([stock_2 stock_3 stock_4 stock_5]);
%     val = c'*inv(R)*c;
%     if(val > 0.6 && val ~= Inf)
%         fprintf('%i, %i, %i, %i, %i, %f\n', ids(i), ids(j), ids(k), ids(l), ids(m), val);
%     end
%     
% end
% 
% pause;

% Multi parallel very instable

%Financials
for sector = sectors(4)' %REMOVE THIS SHIT
    ids = strmatch(sector, stocks.sector);

    %useful to avoid going from 1:end but start from the middle
    ids_length = length(ids);
    for i = 1:ids_length
        fprintf('i = %i\n', i);
        stock_1 = GetPrice(stocks, ids(i));
        
        if(isempty(stock_1))
           continue; 
        end
        
        %[25 26 27]
        %[26 25 27] should be here
        %[26 27 25] we don't care
        %[26 26 25] impossible. move to the next one
        %So i should go through the range and j too with i diff from j
        %So we cheat :) - who cares
        for j = (i+1):ids_length 
            fprintf('j = %i\n', j);
            stock_2 = GetPrice(stocks, ids(j));
            
            if(isempty(stock_2))
                continue; 
            end
            
            for k = (j+1):ids_length
                stock_3 = GetPrice(stocks, ids(k));
                
                if(isempty(stock_3))
                    continue; 
                end
                
                for l = (k+1):ids_length
                    stock_4 = GetPrice(stocks, ids(l));
                    
                    if(isempty(stock_4))
                        continue; 
                    end
                    
                    for m = (l+1):ids_length
                        stock_5 = GetPrice(stocks, ids(m));
                        
                        if(isempty(stock_5))
                            continue; 
                        end
                        
                        cor_mat = corr([stock_1 stock_2 stock_3 stock_4 stock_5]);
                        
%                         if(any(size(cor_mat) ~= [5 5]))
%                            continue; 
%                            %that means one vector contains only NaN
%                         end
                        
                        c = cor_mat(2:5,1);
                        R = corr([stock_2 stock_3 stock_4 stock_5]);
                        val = c'*inv(R)*c;
                        if(val > 0.6 && val ~= Inf)
                            fprintf('%i, %i, %i, %i, %i, %f\n', ids(i), ids(j), ids(k), ids(l), ids(m), val);
                        end
                        
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

