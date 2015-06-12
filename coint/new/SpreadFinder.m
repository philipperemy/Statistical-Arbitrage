
%Works ONLY for lags = 0
function [ spreads ] = SpreadFinder( pp, filename, corr_thres )

    tic;
    stocks       = pp.px;
    stocks_count = size(stocks, 2);
    days_count   = size(stocks, 1);

    f = fopen(filename, 'w');
    s = 'Name 1, Id 1, Sector 1, Name 2, Id 2, Sector 2, Name 3, Id 3, Sector3, P Value No Cointegration, p val r1, pval r2, Corr 12, Corr 23, Corr 13, Corr Ret 12, Corr Ret 23, Corr ret 13, sum h\n';
    fprintf(f, s);
    fprintf(s);

    %init for performance
    m_support = zeros(days_count, stocks_count);
    for i = 1:stocks_count
        m_support(:,i) = isnan(stocks(:,i));
    end

    spreads_count = 0;
    for i = 236:240 %CHANGE THIS SHIT
        for j = (i+1):stocks_count
            for k = (j+1):stocks_count

                if(isequal(m_support(:,k), m_support(:,j)) && isequal(m_support(:,k), m_support(:,i)))

                    stock_1 = GetPrice(pp, i);
                    stock_2 = GetPrice(pp, j);
                    stock_3 = GetPrice(pp, k);

                    try
                        if(isequal(stock_1, stock_2) || isequal(stock_1, stock_3) || isequal(stock_2, stock_3))
                            %move on if at least two stocks are equal
                            continue;
                        end

                        if(corr(stock_1, stock_2) > corr_thres && corr(stock_2, stock_3) > corr_thres && corr(stock_1, stock_3) > corr_thres)
                            
                            [ h, spread, res ] = SpreadConstructor( [i j k], [stock_1, stock_2, stock_3] );
                            
                            if(h)
                                spreads_count = spreads_count + 1;
                                spreads(spreads_count) = spread;

                                ret_1 = diff(log(stock_1));
                                ret_2 = diff(log(stock_2));
                                ret_3 = diff(log(stock_3));

                                s = sprintf('%s, %i, %s, %s, %i, %s, %s, %i, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', ...
                                char(pp.names(i)), i, char(pp.sector(i)),...
                                char(pp.names(j)), j, char(pp.sector(j)),...
                                char(pp.names(k)), k, char(pp.sector(k)),...
                                res.pval_r0_jci, res.pval_r1_jci, res.pval_r2_jci, ...
                                corr(stock_1, stock_2), corr(stock_2, stock_3), corr(stock_1, stock_3), ...
                                corr(ret_1, ret_2), corr(ret_2, ret_3), corr(ret_1, ret_3));
                            
                                fprintf(f, s);
                                fprintf(s);
                            
                            end
                        end
                    catch ex
                        rethrow(ex);
                    end
                end
            end
        end
    end
    toc;
    
    fclose(f);
end

