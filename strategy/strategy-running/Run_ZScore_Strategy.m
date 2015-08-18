addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('..');
load ../../data/spx.mat; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD A SPREAD MODULE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spread_count = length(spreads);
mat 		 = zeros(spread_count , 5 );
zscore_conf  = struct('sell_open', 2, 'sell_close', 0.75, 'buy_open', -2, 'buy_close', -0.5);

for i = 1:spread_count

	Spread 		= spreads(i);
	T			= length(spreads(i));
	pl 			= SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0 );
	cumsum_pl 	= cumsum(pl);
	profit 		= cumsum_pl(end);
	vol 		= std(pl);
	sharpe 		= profit / vol;
	real_sharpe = (mean(pl)/std(pl))*sqrt(252);
	trades 		= sum(pl ~= 0);

	mat(i,1) 	= profit;
	mat(i,2) 	= vol;
	mat(i,3) 	= sharpe;
	mat(i,4) 	= real_sharpe;
	mat(i,5) 	= trades;
	
end