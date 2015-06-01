clear;
load spx.mat;
tuples(1) = 713;
tuples(2) = 949;
tuples(3) = 1187;
Spread = SpreadConstructor(pp, tuples);


%pl = SimpleTradingStrategy( Spread, 1, 500 );
%plot(cumsum(pl));

nstd_range = 0.5:0.1:4.0;
wts_range = 0:1;
wsize_range = 5:1:60;

mat = zeros(length(nstd_range)*length(wts_range)*length(wsize_range),8);
fprintf('ntsd, wts, wsize, profit, vol, sharpe, real_sharpe, trades \n');
i = 1;
for nstd = nstd_range
   for wts = wts_range
      for wsize = wsize_range
   
        boll_conf = struct('wsize', wsize, 'wts', wts, 'nstd', nstd);
        pl = SimpleTradingStrategy( Spread, 1, 6000, boll_conf );
        cumsum_pl = cumsum(pl);
        profit = cumsum_pl(end);
        vol = std(pl);
        sharpe = profit / vol;
        real_sharpe = (mean(pl)/std(pl))*sqrt(252); % 0.7776 http://www.nuclearphynance.com/Show%20Post.aspx?PostIDKey=167056
        trades = sum(pl ~= 0);
        
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

mat = mat(all(~isnan(mat),2),:); % for nan - rows
sorted_mat = sortrows(mat, 6);
l = size(sorted_mat,1);
sorted_mat = sorted_mat((l - 50):l,:);
mean(sorted_mat, 1)

sorted_mat = sorted_mat(all(~isnan(sorted_mat),2),:); % for nan - rows
sorted_mat = sorted_mat(:,all(~isnan(sorted_mat)));   % for nan - columns

%1.2000   60.0000    1.0000   78.8323    0.2682  293.9238  138.0000
boll_conf = struct('wsize', 60, 'wts', 1, 'nstd', 1.2);
pl = SimpleTradingStrategy( Spread, 1, 6000, boll_conf );