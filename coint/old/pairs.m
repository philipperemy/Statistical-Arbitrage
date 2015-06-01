clear;
diary log.txt;
load spx.mat;

stocks = pp.px;
stocks_count = size(stocks, 2);

%pp.px = pp.px(2000:3000,:);

pairs_count = 0;
for i = 1:stocks_count
    for j = 1:stocks_count
        if i ~= j
            
            stock_1 = stocks(:,i);
            nan_stock_1 = isnan(stock_1);
            
            stock_2 = stocks(:,j);
            nan_stock_2 = isnan(stock_2);
            
            if(isnan(stock_1) == isnan(stock_2))
                stock_1(isnan(stock_1)) = [];
                stock_2(isnan(stock_2)) = [];
                
                %pval = Cointegration(stock_1, stock_2);
                try 
                    [~,pval,~,~,~] = jcitest([stock_1, stock_2]);
                    pval = pval.r0;
                    if(pval < 0.05)
                        pairs_count = pairs_count + 1;
                        corr_ret = corr(diff(log(stock_1)), diff(log(stock_2)));
                        corr_prices = corr(stock_1, stock_2);
                        
                        vector_i(pairs_count) = i;
                        vector_j(pairs_count) = j;
                        
                        fprintf('[%i] Id1 = %i, Id2 = %i coint pval = %f, corr ret = %f, corr prices = %f\n', ...
                        pairs_count, i, j, pval, corr_ret, corr_prices);
                    end
                catch
                    warning('Problem.');
                end
            end
        end
    end
end

