function [ portfolio_cumsum, perf_assessment ] = Run_Strategy( spreads, spx, strat_id )

    addpath('../../helpers/');
    addpath('../../coint/deepsearch');
    addpath('../../coint/impl');
    addpath('../../pmcmc');
    addpath('../../filters');
    addpath('../../models');
    addpath('..');
    load ../../data/spx.mat;
    format longg;

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

    initial_bet  = 10000;
    spread_count = length(spreads);
    mat_tr 		 = zeros(spread_count , 5 );
    boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
    zscore_conf  = struct('sell_open', 1, 'sell_close', -1, 'buy_open', -1, 'buy_close', 1);

    %train
    c = 1;
    for i = 1:spread_count
        if(complex_strat)
            try
                fprintf('TRAIN Bollinger Complex, i = %d\n', i);
                s = sprintf('pmcmc_%d.mat', i);
                load(s);
                mat = complex.mat;
                mat_tr(i,1) = mat(1);
                mat_tr(i,2) = mat(2);
                mat_tr(i,3) = mat(3);
                mat_tr(i,4) = mat(4);
                mat_tr(i,5) = mat(5);
                
                complexx(c) = complex;
                c = c + 1;
            catch
            end
        else
            %simple strategy
            Spread 		= spreads(i);
            T           = ceil((1/3)* length(Spread.px));

            if(zscore_strat)
                fprintf('TRAIN ZScore, i = %d\n', i);
                [pl, cum_pro, trds]= SimpleTradingStrategyZScore( pp, Spread, 1, T, zscore_conf, 0, 1, @Strategy_Simulator, initial_bet );
            else
                fprintf('TRAIN Bollinger Simple, i = %d\n', i);
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
    %Select 40 and drop some
    sorted_mat_tr = Rank_Results(mat_tr, 3, 40);
    spreads_ids = flip(sorted_mat_tr(:,5)');
    selected_spreads_ids = [];
    list = {};
    for i = spreads_ids

       if(length(selected_spreads_ids) == 20)
           break;
       end
        
       spread = spreads(i);
       s = sprintf('%s%s%s', char(CleanName(pp.tickers(spread.tuple(1)))), char(CleanName(pp.tickers(spread.tuple(2)))), char(CleanName(pp.tickers(spread.tuple(3)))));
       disp(s);
       s = cellstr(s);
       if(isempty(strmatch(s, list)))
           selected_spreads_ids(end+1) = i;
       end
       
       list(end+1) = s;
    end
    
    spreads_ids      = selected_spreads_ids;
    
    mat_te 		 	 = zeros(length(spreads_ids) , 5);
    portfolio_cumsum = zeros(1);
    c 			 	 = 1;
    for i = spreads_ids

        Spread 		= spreads(i);
        T           = ceil((1/3)* length(Spread.px));

        if(complex_strat)
            for j = 1:length(complexx)
               if(complexx(j).mat(5) == i)
                   spread_vol = complexx(j).model.vol_two_factors_leverage;
                   break;
               end
            end
        else
            spread_vol = 0;
        end

        if(zscore_strat)
            fprintf('TEST ZScore i = %d\n', i);
            [pl, cum_pro, trds]= SimpleTradingStrategyZScore( pp, Spread, T, length(Spread.px), zscore_conf, 0, 1, @Strategy_Simulator, initial_bet );
        else
            if(spread_vol == 0)
                fprintf('TEST Bollinger Simple i = %d\n', i);
            else
                fprintf('TEST Bollinger Complex i = %d\n', i);
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
    figure;
    perf_assessment = PerformanceAssessment(portfolio_cumsum(1:end-2), spx(T:end), initial_bet);
    drawnow;

    DISP_Selected_Triples( pp, spreads, spreads_ids, mat_tr, mat_te );
end