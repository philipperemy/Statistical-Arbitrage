function [ portfolio_cumsum ] = Run_Strategy( spreads, spx, strat_id )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% LOAD A SPREAD MODULE %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ZSCORE_ID           = 3;
    BOLLINGER_CPX       = 2;
    BOLLINGER_SIMPLE    = 1;
    
    switch strat_id
        case ZSCORE_ID
            zscore_strat  = true;
            complex_strat = false;
        case BOLLINGER_CPX
            zscore_strat  = false;
            complex_strat = true;
        case BOLLINGER_SIMPLE
            zscore_strat  = false;
            complex_strat = false;
    end

    if(zscore_strat && complex_strat)
        disp('problem occurred');
        return;
    end

    % load('spreads_7306_8036.mat');
    % load('spx_7306_8036.mat');
    % 
    % load('spreads_8036_8767.mat');
    % load('spx_8036_8767.mat');

    addpath('../../helpers/');
    addpath('../../coint/deepsearch');
    addpath('../../coint/impl');
    addpath('../../pmcmc');
    addpath('../../filters');
    addpath('../../models');
    addpath('..');
    load ../../data/spx.mat;
    format longg;

    initial_bet  = 10000;
    spread_count = length(spreads);
    mat_tr 		 = zeros(spread_count , 5 );
    boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
    zscore_conf  = struct('sell_open', 1, 'sell_close', -1, 'buy_open', -1, 'buy_close', 1);

    %train
    for i = 1:spread_count
        if(complex_strat)
            try
                %fprintf('TRAIN Bollinger Complex, i = %d\n', i);
                s = sprintf('pmcmc_%d.mat', i);
                load(s);
                mat_tr(i,1) = mat(1);
                mat_tr(i,2) = mat(2);
                mat_tr(i,3) = mat(3);
                mat_tr(i,4) = mat(4);
                mat_tr(i,5) = mat(5);
            catch
            end
        else
            %simple strategy

            Spread 		= spreads(i);
            T           = ceil((1/3)* length(Spread.px));

            if(zscore_strat)
                %fprintf('TRAIN ZScore, i = %d\n', i);
                [pl, cum_pro, trds]= SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0, 1, @Strategy_Simulator, initial_bet );
            else
                %fprintf('TRAIN Bollinger Simple, i = %d\n', i);
                spread_vol = 0;
                [pl, cum_pro, trds]= SimpleTradingStrategy( pp, Spread, 1, T, boll_conf, spread_vol, 0, 1, @Strategy_Simulator, initial_bet );
            end

            mat_tr(i,1) = cum_pro(end);
            mat_tr(i,2) = std(pl);
            mat_tr(i,3) = (mean(pl)/std(pl))*sqrt(252);
            mat_tr(i,4) = trds;
            mat_tr(i,5) = i;
        end
    end

    %Rank according to Sharpe Ratio
    sorted_mat_tr = Rank_Results(mat_tr, 3, 20);
    spreads_ids = sorted_mat_tr(:,5)';

    mat_te 		 	 = zeros(length(spreads_ids) , 5 );
    portfolio_cumsum = zeros(1);
    c 			 	 = 1;
    for i = spreads_ids

        Spread 		= spreads(i);
        T           = ceil((1/3)* length(Spread.px));

        if(complex_strat)
            [ pmcmc ]  = Sto_Volatility_Estimator( Spread.px, 65, 600 );
            spread_vol = pmcmc.vol_two_factors_leverage;
        else
            spread_vol = 0;
        end

        if(zscore_strat)
            %fprintf('TEST ZScore i = %d\n', i);
            [pl, cum_pro, trds]= SimpleTradingStrategyZScore( pp, Spread, T, length(Spread.px), zscore_conf, 0, 1, @Strategy_Simulator, initial_bet );
        else
            if(spread_vol == 0)
                %fprintf('TEST Bollinger Simple i = %d\n', i);
            else
                %fprintf('TEST Bollinger Complex i = %d\n', i);
            end

            [pl, cum_pro, trds]= SimpleTradingStrategy( pp, Spread, T, length(Spread.px), boll_conf, spread_vol, 0, 1, @Strategy_Simulator, initial_bet );
        end

        mat_te(c,1)  = cum_pro(end);
        mat_te(c,2)  = std(pl);
        mat_te(c,3)  = (mean(pl)/std(pl))*sqrt(252);
        mat_te(c,4)  = trds;
        mat_te(c,5)  = i;
        try
            mat_te(c,6)  = maxdrawdown(cum_pro);
        catch
            mat_te(c,6) = 1; %lost everything
        end
        
        if(length(portfolio_cumsum) == 1)
           cut  = length(cum_pro);
        else
            cut = length(portfolio_cumsum);
        end
        try
            portfolio_cumsum = portfolio_cumsum + cum_pro(1:cut);
        catch
            cum_pro(end+1) = cum_pro(end);
            portfolio_cumsum = portfolio_cumsum + cum_pro(1:cut);
        end
        c = c + 1;
    end
    portfolio_cumsum = portfolio_cumsum / length(spreads_ids);
    PerformanceAssessment(portfolio_cumsum, spx(T:end), initial_bet);
    drawnow;

    DISP_Selected_Triples( pp, spreads, spreads_ids, mat_tr, mat_te );
end