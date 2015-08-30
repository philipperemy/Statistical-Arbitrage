clear;clc;
addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('../../pmcmc');
addpath('../../filters');
addpath('../../models');
addpath('..');
load ../../data/spx.mat;
format longg;

load('spreads_5114_5845.mat');
load('spx_5114_5845.mat');

mat_a = zeros(length(spreads),5);
for i = 1:length(spreads)
    try
        s = sprintf('pmcmc_%d.mat', i);
        load(s);
        mat = complex.mat;
        mat_a(i,1) = mat(1);
        mat_a(i,2) = mat(2);
        mat_a(i,3) = mat(3);
        mat_a(i,4) = mat(4);
        mat_a(i,5) = mat(5);
    catch
        
    end
end

sorted_mat_a = Rank_Results(mat_a, 3, 20);
spreads_ids = sorted_mat_a(:,5)';
%sr = mat_a(:,3); sr(isnan(sr)) = []; sr(sr == 0) = []; mean(sr)
mat_t 		 = zeros(length(spreads_ids) , 5 );
fprintf('testing\n');
c = 1;
portfolio_cumsum = zeros(1);
boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
initial_bet  = 10000;

for i = spreads_ids
    fprintf('i = %d\n', i);
	Spread 		= spreads(i);
    T           = ceil((1/3)* length(Spread.px));
    [ pmcmc ]   = Sto_Volatility_Estimator( Spread.px, 65, 3000 );
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
PerformanceAssessment(portfolio_cumsum, spx(T:end), initial_bet)
portfolio_cumsum_sv = portfolio_cumsum;


%double dataset not that goood.
tuples = [];
for i = spreads_ids
   spread = spreads(i);
   
   c = find(mat_t(:,5)==i);
   sro = mat_t(c,3);
   tr = mat_t(c,4);
   mdd = mat_t(c,6);
   net = mat_t(c,1);

   sri = mat_a(i,3);
   s = sprintf('%d, %s & %s & %s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %i \\\\', i, char(CleanName(pp.tickers(spread.tuple(1)))), char(CleanName(pp.tickers(spread.tuple(2)))), char(CleanName(pp.tickers(spread.tuple(3)))),...
       spread.beta(1), spread.beta(2), spread.beta(3), sri, sro, net, mdd, tr);
   tuples(end+1) = spread.tuple(1);
   tuples(end+1) = spread.tuple(2);
   tuples(end+1) = spread.tuple(3);
   disp(s);
end

tuples = unique(tuples);
for i = 1:length(pp.tickers)
   if(sum(find(i==tuples)~=0))
       s = sprintf('%s & %s & %s \\\\', char(CleanName(pp.tickers(i))), char(pp.names(i)), char(pp.sector(i)));
       disp(s);
   end
end
