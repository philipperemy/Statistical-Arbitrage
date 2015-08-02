run('Prepare_Stocks.m');

vec = zeros(1,100);

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
                    if(val < 0.999)
                        vec = Update_Map(vec, val);
                        if(val > max_R) %avoid 100% correlation
                            max_R = val;
                            fprintf('%i (%i), %i (%i), %i (%i), %i (%i), %f\n', i, ids(i), j, ids(j), k, ids(k), l, ids(l), max_R);
                        end
                    end
                end
            end
        end
    end
end

% compute R^2 on stationary series
%155, 157, 309, 341, 0.999931
% 5, 6, 20, 21
% stock_1 = stocks.ret(:,5);
% stock_2 = stocks.ret(:,6);
% stock_3 = stocks.ret(:,20);
% stock_4 = stocks.ret(:,21);
% mdl = fitlm([stock_2 stock_3 stock_4], stock_1, 'Intercept', false);
% val = mdl.Rsquared.Ordinary;