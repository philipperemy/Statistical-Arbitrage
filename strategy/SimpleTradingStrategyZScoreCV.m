addpath('../helpers/');
load ../data/spx.mat; 
spx = csvread('../data/spx.csv');
addpath('../coint/deepsearch');
addpath('../coint/impl');

ids = [281 289 726];
[ ~, Spread, ~ ] = SpreadConstructorLog(ids, log(GetPriceArray(pp, ids)));
Spread.px = exp(Spread.px); %exp price

zscore_conf = struct('sell_open', 2, 'sell_close', 0.75, 'buy_open', -2, 'buy_close', -0.5);
ZScore(Spread.px, zscore_conf);


sell_open_range   = 1.5:0.1:2.5;
sell_close_range  = 0.5:0.1:0.8;
buy_open_range    = -2.5:0.1:-1.5;
buy_close_range   = -0.7:0.1:-0.5;

mat = zeros(length(sell_open_range)*length(sell_close_range)*length(buy_open_range)*length(buy_close_range),9);
fprintf('sell_open, sell_close, buy_open, buy_close, profit, vol, sharpe, real_sharpe, trades \n');
i = 1;
T = length(Spread.px);
for sell_open = sell_open_range
   for sell_close = sell_close_range
      for buy_open = buy_open_range
        for buy_close = buy_close_range
            
            zscore_conf = struct('sell_open', sell_open, 'sell_close', sell_close, 'buy_open', buy_open, 'buy_close', buy_close);
            pl = SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0 );
            cumsum_pl = cumsum(pl);
            profit = cumsum_pl(end);
            vol = std(pl);
            sharpe = profit / vol;
            real_sharpe = (mean(pl)/std(pl))*sqrt(252); % 0.7776 http://www.nuclearphynance.com/Show%20Post.aspx?PostIDKey=167056
            trades = sum(pl ~= 0);

            fprintf('%f, %f, %f, %f, %f, %f, %f, %f, %i\n', sell_open, sell_close, buy_open, buy_close, profit, vol, sharpe, real_sharpe, trades);

            mat(i,1) = sell_open;
            mat(i,2) = sell_close;
            mat(i,3) = buy_open;
            mat(i,4) = buy_close;
            mat(i,5) = profit;
            mat(i,6) = vol;
            mat(i,7) = sharpe;
            mat(i,8) = real_sharpe;
            mat(i,9) = trades;

            i = i + 1;
        end
      end
   end
end

mat = mat(all(~isnan(mat),2),:); % for nan - rows
sorted_mat = sortrows(mat, 5);
l = size(sorted_mat,1);
sorted_mat = sorted_mat((l - 50):l,:);
mean(sorted_mat, 1)

sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns

zscore_conf = struct('sell_open', 1.8, 'sell_close', 0.6, 'buy_open', -1.7, 'buy_close', -0.5);
pl = SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0 );
cumsum_pl = cumsum(pl);