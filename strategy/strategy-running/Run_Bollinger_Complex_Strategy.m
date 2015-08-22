addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('../../pmcmc');
addpath('../../filters');
addpath('../../models');
addpath('..');
load ../../data/spx.mat;
format longg;

load('spreads_6575_7306.mat');
load('spx_6575_7306.mat');

%%%%%%%%%%%%%%
% Run Simple before and load spreads_ids and spreads
%%%%%%%%%%%%%%

% N = length(spreads(1).px)+100;
% for spr_id = spreads_ids
%     spread      = spreads(spr_id);
%     [ pmcmc ]   = Sto_Volatility_Estimator( spread.px, 100, 3000 );
%     fname       = sprintf('pmcmc_%d.mat', spr_id); 
%     save(fname, 'pmcmc');
% end



mat_t 		 = zeros(length(spreads_ids) , 5 );
fprintf('testing\n');
c = 1;
portfolio_cumsum = zeros(1);
for i = spreads_ids
    fprintf('i = %d\n', i);
	Spread 		= spreads(i);
    T           = ceil((1/3)* length(Spread.px));
    [ pmcmc ]   = Sto_Volatility_Estimator( Spread.px, 65, 6000 );
	spread_vol = pmcmc.vol_two_factors_leverage;
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

%Sometimes you can just remove the end if it's too bad.
PerformanceAssessment(portfolio_cumsum, spx(T:end), initial_bet)