run('Init.m');
stock_ids = [713 949 1187];
stocks_prices = [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ];
[ h, Spread, ~ ] = SpreadConstructor(stock_ids, [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ]);
% 
% nstd_range = 1.0:0.1:3.0;
% wts_range = 0:1;
% wsize_range = 5:1:100;
% 
% mat = zeros(length(nstd_range)*length(wts_range)*length(wsize_range),8);
% fprintf('ntsd, wts, wsize, profit, vol, sharpe, real_sharpe, trades \n');
% i = 1;
% for nstd = nstd_range
%    for wts = wts_range
%       for wsize = wsize_range
%    
%         boll_conf = struct('wsize', wsize, 'wts', wts, 'nstd', nstd);
%         pl = SimpleTradingStrategy( Spread, 1, 6000, boll_conf );
%         cumsum_pl = cumsum(pl);
%         profit = cumsum_pl(end);
%         vol = std(pl);
%         sharpe = profit / vol;
%         real_sharpe = (mean(pl)/std(pl))*sqrt(252); % 0.7776 http://www.nuclearphynance.com/Show%20Post.aspx?PostIDKey=167056
%         trades = sum(pl ~= 0);
%         
%         fprintf('%f, %i, %i, %f, %f, %f, %f, %i\n', nstd, wts, wsize, profit, vol, sharpe, real_sharpe, trades);
%         
%         mat(i,1) = nstd;
%         mat(i,2) = wsize;
%         mat(i,3) = wts;
%         mat(i,4) = profit;
%         mat(i,5) = vol;
%         mat(i,6) = sharpe;
%         mat(i,7) = real_sharpe;
%         mat(i,8) = trades;
%         
%         i = i + 1;
%       end
%    end
% end
% 
% mat = mat(all(~isnan(mat),2),:); % for nan - rows
% sorted_mat = sortrows(mat, 7);
% l = size(sorted_mat,1);
% sorted_mat = sorted_mat((l - 50):l,:);
% mean(sorted_mat, 1)
% 
% sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
% sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns

boll_conf = struct('wsize', 74, 'wts', 1, 'nstd', 2.5);
[pl,balance_cum] = SimpleTradingStrategy( pp, Spread, 1, 6100, boll_conf );
plot(balance_cum);

PerformanceAssessment(balance_cum, spx);
