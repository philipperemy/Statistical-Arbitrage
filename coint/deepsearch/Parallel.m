run('../Init.m');
stock_ids = [713 949 1187];
stocks_prices = [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ];
[ h, Spread, ~ ] = SpreadConstructor(stock_ids, [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ]);

nstd_range = 1.0:0.1:3.0;
wts_range = 0:1;
wsize_range = 5:1:100;

mat = zeros(length(nstd_range)*length(wts_range)*length(wsize_range),8);
fprintf('ntsd, wts, wsize, profit, vol, sharpe, real_sharpe, trades \n');

tic;
I = combvec(nstd_range, wts_range, wsize_range)';
O = zeros(1, length(I));
parfor j = 1:length(I);
    wsize = I(j,3);
    wts = I(j,2);
    nstd = I(j,1);
    boll_conf = struct('wsize', wsize, 'wts', wts, 'nstd', nstd);
    pl = SimpleTradingStrategy( pp, Spread, 1, 6000, boll_conf, 0 );
    cumsum_pl = cumsum(pl);
    O(j) = cumsum_pl(end);
end
toc;
fprintf('end\n');
pause;

[O1, O2, O3] = CrossValidationParallel( pp, Spread, @SimpleTradingStrategy, nstd_range, wts_range, wsize_range );
fprintf('end 2\n');
pause;

tic;
i = 1;
for nstd = nstd_range
   for wts = wts_range
      for wsize = wsize_range

        boll_conf = struct('wsize', wsize, 'wts', wts, 'nstd', nstd);
        pl = SimpleTradingStrategy( pp, Spread, 1, 6000, boll_conf, 0 );
        cumsum_pl = cumsum(pl);
        profit = cumsum_pl(end);
        vol = std(pl);
        sharpe = profit / vol;
        real_sharpe = (mean(pl)/std(pl))*sqrt(252); % 0.7776 http://www.nuclearphynance.com/Show%20Post.aspx?PostIDKey=167056
        trades = sum(pl ~= 0);
        fprintf('length pl %i\n', length(pl));
        fprintf('%f, %i, %i, %f, %f, %f, %f, %i\n', nstd, wts, wsize, profit, vol, sharpe, real_sharpe, trades);
        
        mat(i,1) = nstd;
        mat(i,2) = wsize;
        mat(i,3) = wts;
        mat(i,4) = profit;
        mat(i,5) = vol;
        mat(i,6) = sharpe;
        mat(i,7) = real_sharpe;
        mat(i,8) = trades;
        
        i = i + 1;
      end
   end
end
toc;

mat = mat(all(~isnan(mat),2),:); % for nan - rows
sorted_mat = sortrows(mat, 7);
l = size(sorted_mat,1);
sorted_mat = sorted_mat((l - 50):l,:);
mean(sorted_mat, 1)

sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns

boll_conf = struct('wsize', 74, 'wts', 1, 'nstd', 2.5);
[pl,balance_cum] = SimpleTradingStrategy( pp, Spread, 1, 6100, boll_conf );
plot(balance_cum);

PerformanceAssessment(balance_cum, spx);


%portfolio
[~,b1] = SimpleTradingStrategy( pp, Spread, 1, 6000, struct('wsize', sorted_mat(end,2), 'wts', sorted_mat(end,3), 'nstd', sorted_mat(end,1)) );
[~,b2] = SimpleTradingStrategy( pp, Spread, 1, 6000, struct('wsize', sorted_mat(end-1,2), 'wts', sorted_mat(end-1,3), 'nstd', sorted_mat(end-1,1)) );
[~,b3] = SimpleTradingStrategy( pp, Spread, 1, 6000, struct('wsize', sorted_mat(end-2,2), 'wts', sorted_mat(end-2,3), 'nstd', sorted_mat(end-2,1)) );
[~,b4] = SimpleTradingStrategy( pp, Spread, 1, 6000, struct('wsize', sorted_mat(end-3,2), 'wts', sorted_mat(end-3,3), 'nstd', sorted_mat(end-3,1)) );
[~,b5] = SimpleTradingStrategy( pp, Spread, 1, 6000, struct('wsize', sorted_mat(end-4,2), 'wts', sorted_mat(end-4,3), 'nstd', sorted_mat(end-4,1)) );

plot([b1' b2' b3' b4' b5']);

y = [b1' b2' b3' b4' b5'];
w = repmat(0.2, 1, 5);
portfolio = y*w';
pl = diff(portfolio);
real_sharpe = (mean(pl)/std(pl))*sqrt(252);
real_sharpe

%http://katrinaeg.com/simulated-annealing.html