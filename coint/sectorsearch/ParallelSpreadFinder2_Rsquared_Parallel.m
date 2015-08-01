
run('../Init.m');
addpath('../../helpers');
diary off;
delete 'Financials.txt';
diary 'Financials.txt';
tic;

twobytwo            = 1;
year_start_dates    = Extract_Years(pp, twobytwo);
last_start          = year_start_dates(end-1);
last_end            = year_start_dates(end);

[pp_train_2013, pp_test_2013] = Create_Sets( pp, last_start, last_end, 0.7 );
clear pp; %just to avoid any bugs

stocks          =   Prepare_Set( pp_train_2013 );
sectors         =   unique(stocks.sector);
stocks          =   rmfield(stocks, 'px_clean');
stocks          =   rmfield(stocks, 'w');

% compute R^2 on stationary series

%[25 26 27]
%[26 25 27] should be here
%[26 27 25] we don't care
%[26 26 25] impossible. move to the next one
%So i should go through the range and j too with i diff from j
%So we cheat :) - who cares

% max_R = 0;
% %Consumer Staples
% for sector = sectors(2)' %REMOVE THIS SHIT
%     ids = strmatch(sector, stocks.sector);
%     stock_etf = diff(mean(normc(GetPriceArray(stocks, ids)),2));
%     
%     %useful to avoid going from 1:end but start from the middle
%     ids_length = length(ids);
%     for i = 1:ids_length
%         fprintf('i = %i\n', i);
%         stock_1 = stocks.ret(:,i);
% 
%         for j = (i+1):ids_length 
%             fprintf('j = %i\n', j);
%             stock_2 = stocks.ret(:,j);
%             
%             for k = (j+1):ids_length
%                 stock_3 = stocks.ret(:,k);
%                 
%                 for l = (k+1):ids_length
%                     stock_4 = stocks.ret(:,l);
%    
%                     A = [stock_etf stock_1 stock_2 stock_3 stock_4];
%                     cor_mat = corr(A);
%                     c = cor_mat(2:5,1);
%                     A(:,1) = [];
%                     R = corr(A);
%                     val = c'*inv(R)*c;
%                     if(val > max_R)
%                         max_R = val;
%                         fprintf('ETF, %i, %i, %i, %i, %f\n', ids(i), ids(j), ids(k), ids(l), max_R);
%                     end
%                 end
%             end
%         end
%     end
% end


%%%WITHOUT ETF

%155, 157, 309, 341, 0.999931
% 5, 6, 20, 21
% stock_1 = stocks.ret(:,5);
% stock_2 = stocks.ret(:,6);
% stock_3 = stocks.ret(:,20);
% stock_4 = stocks.ret(:,21);
% mdl = fitlm([stock_2 stock_3 stock_4], stock_1, 'Intercept', false);
% val = mdl.Rsquared.Ordinary;


%Financials
for sector = sectors'
    ids = strmatch(sector, stocks.sector);
    max_R = 0;
    fprintf('SECTOR:');
    disp(sector);
    fprintf('Index 1, Index 2, Index 3, Index 4, R^2');
    ids_length = length(ids);
    for i = 1:ids_length
        stock_1 = stocks.ret(:,i);

        for j = (i+1):ids_length 
            stock_2 = stocks.ret(:,j);
            
            for k = (j+1):ids_length
                stock_3 = stocks.ret(:,k);
                
                for l = (k+1):ids_length
                    stock_4 = stocks.ret(:,l);
                    
                    val = Compute_R_2( stock_1, stock_2, stock_3, stock_4 );  
                    if(val > max_R && val < 0.999) %avoid 100% correlation
                        max_R = val;
                        fprintf('%i (%i), %i (%i), %i (%i), %i (%i), %f\n', i, ids(i), j, ids(j), k, ids(k), l, ids(l), max_R);
                    end
                end
            end
        end
    end
end

diary off;
diary on;
