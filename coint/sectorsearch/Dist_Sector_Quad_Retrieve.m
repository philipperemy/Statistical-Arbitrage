function [ quadruples ] = Dist_Sector_Quad_Retrieve( stocks, sector, thres, debug )
    %
    % HIGH PERFORMANCE 1 (i+1) (j+1) (k+1) instead of 1 1 1 1
    %
    quadruples = [];
    stocks_sector = RetrieveStocksSector(stocks, sector);
    stocks_count = length(stocks_sector.sector);
    fprintf('SECTOR:');
    disp(sector);
    fprintf('Index 1, Index 2, Index 3, Index 4, R^2\n');
    fprintf('THRESHOLD:');
    disp(thres);
    total_iterations = stocks_count^4;
    fprintf('expected number of iterations: %i\n', total_iterations);
    count = 1;
    for i = 1:stocks_count
        stock_1 = stocks_sector.ret(:,i);

        for j = (i+1):stocks_count 
            stock_2 = stocks_sector.ret(:,j);

            for k = (j+1):stocks_count
                stock_3 = stocks_sector.ret(:,k);

                for l = (k+1):stocks_count
                    stock_4 = stocks_sector.ret(:,l);

                    if(IdsAreDifferent(i,j,k,l))

                        val = Compute_R_2( stock_1, stock_2, stock_3, stock_4 );
                        if(val < 0.999 && val > thres)
                            fprintf('%i (%s), %i (%s), %i (%s), %i (%s), %f\n', i,  char(stocks_sector.names(i)), ...
                                                                                j, char(stocks_sector.names(j)), ...
                                                                                k, char(stocks_sector.names(k)), ...
                                                                                l, char(stocks_sector.names(l)), ...
                                                                                val);
                            quadruples = [quadruples; [i j k l val]];
                            count = count + 1;
                            if(count == 100 && exist('debug', 'var'))
                                return;
                            end
                        end
                    end
                end
            end
        end
    end

end

