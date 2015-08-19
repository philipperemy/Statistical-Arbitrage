addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('..');
load ../../data/spx.mat;
format longg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD A SPREAD MODULE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spread_count = length(spreads);
mat 		 = zeros(spread_count , 5 );
%zscore_conf  = struct('sell_open', 2, 'sell_close', 0.75, 'buy_open', -2, 'buy_close', -0.5);
zscore_conf  = struct('sell_open', 1, 'sell_close', -1, 'buy_open', -1, 'buy_close', 1);

for i = 1:spread_count
    fprintf('i = %d\n', i);
	Spread 		= spreads(i);
	T			= length(Spread.px);
	[pl, cum_pro]= SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0, 1 );
	vol 		= std(pl);
	sharpe 		= profit / vol;
	real_sharpe = (mean(pl)/std(pl))*sqrt(252);
	trades 		= sum(pl ~= 0);

	mat(i,1) 	= cum_pro(end);
	mat(i,2) 	= vol;
	mat(i,3) 	= real_sharpe;
	mat(i,4) 	= trades;
    mat(i,5)    = i;
	
end

sorted_mat       = Rank_Results( mat, 3, 50);
spreads_ids      = sorted_mat(end-10:end,5)';
portfolio_cumsum = Cumsum_PL( pp, spreads, spreads_ids, zscore_conf );