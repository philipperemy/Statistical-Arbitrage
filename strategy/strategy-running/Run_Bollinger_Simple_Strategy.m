%clear;clc;
load('spreads_6575_7306.mat');
load('spx_6575_7306.mat');

% load('spreads_7306_8036.mat');
% load('spx_7306_8036.mat');
% 
% load('spreads_8036_8767.mat');
% load('spx_8036_8767.mat');

addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('..');
load ../../data/spx.mat;
format longg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD A SPREAD MODULE %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%portfolio = [portfolio portfolio_cumsum];
initial_bet  = 10000;

spread_count = length(spreads);
mat 		 = zeros(spread_count , 5 );
boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
spread_vol   = 0;
%train
for i = 1:spread_count
    fprintf('i = %d\n', i);
    Spread 		= spreads(i);
    T           = ceil((1/3)* length(Spread.px));
    [pl, cum_pro, trds]= SimpleTradingStrategy( pp, Spread, 1, T, boll_conf, spread_vol, 0, 1, @Strategy_Simulator, initial_bet );
    vol          = std(pl);
    sharpe_ratio = (mean(pl)/std(pl))*sqrt(252);
	trades 		 = trds;

	mat(i,1) 	= cum_pro(end);
	mat(i,2) 	= vol;
	mat(i,3) 	= sharpe_ratio;
	mat(i,4) 	= trades;
    mat(i,5)    = i;
end

sorted_mat  = Rank_Results(mat, 3, 50);
spreads_ids = sorted_mat(end-19:end,5)';

mat_t 		 = zeros(length(spreads_ids) , 5 );
fprintf('testing\n');
c = 1;
portfolio_cumsum = zeros(1);
for i = spreads_ids
    fprintf('i = %d\n', i);
	Spread 		= spreads(i);
    T           = ceil((1/3)* length(Spread.px));
	[pl, cum_pro, trds]= SimpleTradingStrategy( pp, Spread, T, length(Spread.px), boll_conf, spread_vol, 0, 1, @Strategy_Simulator, initial_bet );
	vol 		= std(pl);
	sharpe_ratio= (mean(pl)/std(pl))*sqrt(252);

	mat_t(c,1) 	= cum_pro(end);
	mat_t(c,2) 	= vol;
	mat_t(c,3) 	= sharpe_ratio;
	mat_t(c,4) 	= trds;
    mat_t(c,5)  = i;
    mat_t(c,6)  = maxdrawdown(cum_pro);
    portfolio_cumsum = portfolio_cumsum + cum_pro;
    c = c + 1;
end
portfolio_cumsum = portfolio_cumsum / length(spreads_ids);
portfolio_cumsum_b = portfolio_cumsum;
PerformanceAssessment(portfolio_cumsum, spx(T:end), initial_bet)

DISP_Selected_Triples( pp, spreads, spreads_ids, mat, mat_t );