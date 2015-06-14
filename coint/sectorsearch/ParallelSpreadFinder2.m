run('../Init.m');
diary off;
delete 'log.txt';
diary 'log.txt';
tic;

stocks          =   pp.px;
spreads_count   =   1;
sectors         =   unique(pp.sector);

[~,unique_idx]  =   unique(pp.px_clean','rows');
pp.px_clean     =   pp.px_clean(:,unique_idx);
pp.sector       =   pp.sector(unique_idx);    
pp.tickers      =   pp.tickers(unique_idx);
pp.names        =   pp.names(unique_idx);
pp.px           =   pp.px(:,unique_idx);
pp.w            =   pp.w(:,unique_idx);

pp.px_clean(1,:)=   [];

pp.px(1,:)      =   [];
pp.dt(1)        =   [];

%Financials
for sector = sectors(4)' %REMOVE THIS SHIT
    ids = strmatch(sector, pp.sector);
    ids_length = length(ids);
    parfor i = 1:ids_length
        fprintf('i = %i\n', i);
        stock_1 = pp.px_clean(:,ids(i));
        
        for j = (i+1):ids_length
            stock_2 = pp.px_clean(:,ids(j));
            
            for k = (j+1):ids_length
                stock_3 = pp.px_clean(:,ids(k));
                
                for l = (k+1):ids_length
                    stock_4 = pp.px_clean(:,ids(l));
                    
                    for m = (l+1):ids_length
                        stock_5 = pp.px_clean(:,ids(m));
                        
                        cor = corr([stock_1 stock_2 stock_3 stock_4 stock_5]);
                        mean_cor = mean(mean(cor));
                        if(mean_cor > 0.83)
                            
                            ret_1 = diff(log(stock_1));
                            ret_2 = diff(log(stock_2));
                            ret_3 = diff(log(stock_3));
                            ret_4 = diff(log(stock_4));
                            ret_5 = diff(log(stock_5));

                            if(any(isnan(ret_1)) || any(isnan(ret_2)) || any(isnan(ret_3)) || any(isnan(ret_4)) || any(isnan(ret_5)))
                                continue;
                            end

                            corr_ret = corr([ret_1 ret_2 ret_3 ret_4 ret_5]);
                            mean_cor_ret = mean(mean(corr_ret));
                            fprintf('%i, %i, %i, %i, %i, %f, %f\n', i, j, k, l, m, mean_cor, mean_cor_ret);

                        end
                    end
                end
            end
        end
    end
end

diary off;
diary on;
