run('../Init.m');
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

%Consumer Discretionary
for sector = sectors(1)' %REMOVE THIS SHIT
    ids = strmatch(sector, pp.sector);
    ids_length = length(ids);
    for i = 1:ids_length
        fprintf('i = %i\n', i);
        for j = (i+1):ids_length
            for k = (j+1):ids_length
                for l = (k+1):ids_length
                    for m = (l+1):ids_length
                        
                        s1 = ids(i);
                        s2 = ids(j);
                        s3 = ids(k);
                        s4 = ids(l);
                        s5 = ids(m);

                        stock_1 = pp.px_clean(:,s1);
                        stock_2 = pp.px_clean(:,s2);
                        stock_3 = pp.px_clean(:,s3);
                        stock_4 = pp.px_clean(:,s4);
                        stock_5 = pp.px_clean(:,s5);

                        cor = corr([stock_1 stock_2 stock_3 stock_4 stock_5]);
                        mean_cor = mean(mean(cor));
                        if(mean_cor > 0.83)
                            fprintf('%i, %i, %i, %i, %i\n', i, j, k, l, m);
                        end
                    end
                end
            end
        end
    end
end

% stock_ids = [s1 s2 s3 s4 s5];
% [ h, spread, res ] =  SpreadConstructor2( stock_ids, stock_prices );
% 
% if(h)
%     spreads(spreads_count) = spread;
%     fprintf('%i\n', spreads_count);
%     spreads_count = spreads_count + 1;
% end