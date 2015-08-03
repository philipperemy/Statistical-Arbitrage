function [ dist_mat ] = Dist_Sector_Quad( stocks, sectors )
    dist_mat = zeros(length(sectors), 100);
    sectorid = 0;
    for sector = sectors'
        sectorid = sectorid + 1;
        stocks_sector = RetrieveStocksSector(stocks, sector);
        max_R = 0;
        stocks_count = length(stocks_sector.sector);
        counter = 0; tic;
        fprintf('SECTOR:');
        disp(sector);
        fprintf('Index 1, Index 2, Index 3, Index 4, R^2\n');
        total_iterations = stocks_count^4;
        fprintf('expected number of iterations: %i\n', total_iterations);
        for i = 1:stocks_count
            stock_1 = stocks_sector.ret(:,i);

            for j = 1:stocks_count 
                stock_2 = stocks_sector.ret(:,j);

                for k = 1:stocks_count
                    stock_3 = stocks_sector.ret(:,k);

                    for l = 1:stocks_count
                        stock_4 = stocks_sector.ret(:,l);

                        if(IdsAreDifferent(i,j,k,l))

                            counter = counter + 1;
                            if(~rem(counter, 10000))
                                estimated_elapsed_time = toc*total_iterations/10000; tic;
                                fprintf('[debug] counter = %i, estimated time (sec): %i\n', counter, estimated_elapsed_time);
                            end

                            val = Compute_R_2( stock_1, stock_2, stock_3, stock_4 );
                            if(val < 0.999)
                                dist_mat(sectorid, :) = Update_Map(dist_mat(sectorid, :), val);
                                if(val > max_R)
                                    max_R = val;
                                    fprintf('%i (%s), %i (%s), %i (%s), %i (%s), %f\n', i,  char(stocks_sector.names(i)), ...
                                                                                        j, char(stocks_sector.names(j)), ...
                                                                                        k, char(stocks_sector.names(k)), ...
                                                                                        l, char(stocks_sector.names(l)), ...
                                                                                        max_R);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end

