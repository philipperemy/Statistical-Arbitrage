run('Prepare_Stocks.m');

for sector = sectors'
    ids = strmatch(sector, stocks.sector);
    max_R = 0;
    stock_etf = diff(mean(normc(GetPriceArray(stocks, ids)),2));
    fprintf('SECTOR:');
    disp(sector);
    fprintf('ETF, Index 1, Index 2, Index 3, Index 4, R^2');
    ids_length = length(ids);
    for i = 1:ids_length
        stock_1 = stocks.ret(:,i);

        for j = (i+1):ids_length
            stock_2 = stocks.ret(:,j);
            
            for k = (j+1):ids_length
                stock_3 = stocks.ret(:,k);
                
                for l = (k+1):ids_length
                    stock_4 = stocks.ret(:,l);

                    val = Compute_R_2_ETF(stock_etf, stock_1, stock_2, stock_3, stock_4);
                    if(val > max_R && val < 0.999) %avoid 100% correlation
                        max_R = val;
                        fprintf('ETF, %i, %i, %i, %i, %f\n', ids(i), ids(j), ids(k), ids(l), max_R);
                    end
                end
            end
        end
    end
end