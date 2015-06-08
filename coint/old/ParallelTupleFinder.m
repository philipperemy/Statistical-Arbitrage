clear;
diary log.txt;
load spx.mat;

stocks = pp.px;
stocks_count = size(stocks, 2);
days_count = size(stocks, 1);

fprintf('Name 1, Id 1, Sector 1, Name 2, Id 2, Sector 2, Name 3, Id 3, Sector3, P Value No Cointegration, p val r1, pval r2, Corr 12, Corr 23, Corr 13, Corr Ret 12, Corr Ret 23, Corr ret 13\n');

%burn
m = zeros(days_count, stocks_count);
for i = 1:stocks_count
    m(:,i) = isnan(stocks(:,i));
end

tic;
thres_corr = 0.80;
ids = zeros(1,1);

parfor i = 1:stocks_count
    for j = (i+1):stocks_count
        for k = (j+1):stocks_count
 
            if(isequal(m(:,k), m(:,j)) && isequal(m(:,k), m(:,i)))

                stock_1 = stocks(:,i);
                stock_2 = stocks(:,j);
                stock_3 = stocks(:,k);

                stock_1(isnan(stock_1)) = [];
                stock_2(isnan(stock_2)) = [];
                stock_3(isnan(stock_3)) = [];

                try
                    if(isequal(stock_1, stock_2) || isequal(stock_1, stock_3) || isequal(stock_2, stock_3))
                        %continue if at least two stocks are equal
                        continue; 
                    end

                    if(corr(stock_1, stock_2) > thres_corr && corr(stock_2, stock_3) > thres_corr && corr(stock_1, stock_3) > thres_corr)
                        [~,pval,~,~,~] = jcitest([stock_1, stock_2, stock_3]); %Check lags as well.
                        pval_r1 = pval.r1;
                        pval_r2 = pval.r2;
                        pval_r0 = pval.r0;
                        if(pval_r0 < 0.05)
                            
                            ret_1 = diff(log(stock_1));
                            ret_2 = diff(log(stock_2));
                            ret_3 = diff(log(stock_3));

                            fprintf('%s, %i, %s, %s, %i, %s, %s, %i, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', ...
                            char(pp.names(i)), i, char(pp.sector(i)),...
                            char(pp.names(j)), j, char(pp.sector(j)),...
                            char(pp.names(k)), k, char(pp.sector(k)),...
                            pval_r0, pval_r1, pval_r2, ...
                            corr(stock_1, stock_2), corr(stock_2, stock_3), corr(stock_1, stock_3), ...
                            corr(ret_1, ret_2), corr(ret_2, ret_3), corr(ret_1, ret_3));

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