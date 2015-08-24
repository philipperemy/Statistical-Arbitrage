function [ ] = PMCMC_Scripting( seq_i )

    load('spreads_6575_7306.mat');
    load('spx_6575_7306.mat');
    
    addpath('../../helpers/');
    addpath('../../coint/deepsearch');
    addpath('../../coint/impl');
    addpath('../../models');
    addpath('../../pmcmc');
    addpath('../../filters');
    addpath('..');
    load ../../data/spx.mat;
    format longg;

    initial_bet  = 10000;

    spread_count = length(spreads);
    boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
    %train
    for i = seq_i
        fprintf('i = %d\n', i);
        Spread 		= spreads(i);
        T           = ceil((1/3)* length(Spread.px));
        [ pmcmc ]   = Sto_Volatility_Estimator( Spread.px, 65, 1500 );
        spread_vol  = pmcmc.vol_two_factors_leverage;
        [pl, cum_pro, trds]= SimpleTradingStrategy( pp, Spread, 1, T, boll_conf, spread_vol, 0, 1, @Strategy_Simulator, initial_bet );
        vol          = std(pl);
        sharpe_ratio = (mean(pl)/std(pl))*sqrt(252);
        trades 		 = trds;

        mat = zeros(1,5);
        mat(1) 	= cum_pro(end);
        mat(2) 	= vol;
        mat(3) 	= sharpe_ratio;
        mat(4) 	= trades;
        mat(5)  = i;
        fname       = sprintf('pmcmc_%d.mat', i ); 
        save(fname, 'mat');
    end

end

