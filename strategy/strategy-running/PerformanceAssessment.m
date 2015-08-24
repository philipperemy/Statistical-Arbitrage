function [ s ] = PerformanceAssessment( strat_porfolio, spx_index, initial_bet )

    if(~exist('initial_bet', 'var'))
        initial_bet = strat_porfolio(1);
    end
    
    T = length(spx_index)-length(strat_porfolio)+1;
    spx_portfolio = ((initial_bet/spx_index(T))*spx_index(T:end));
    
    if(strat_porfolio(1) ~= initial_bet)
       strat_porfolio =  ((initial_bet/strat_porfolio(1))*strat_porfolio);
    end
    
    pl_spx = diff(spx_portfolio);
    pl_strat = diff(strat_porfolio);

    %if the risk free rate is 2% = *1.02 per year
    risk_free_rate = 1.02;
    daily_percentage = risk_free_rate^(1/365);
    daily_risk_free_rate = repmat(daily_percentage, 1, length(strat_porfolio));
    risk_free_rate_portfolio = initial_bet*cumprod(daily_risk_free_rate)';
    pl_risk_free = diff(risk_free_rate_portfolio);
   
    pl_strat_risk_free = pl_strat - mean(pl_risk_free);
    sharpe_ratio_free = (mean(pl_strat_risk_free) / std(pl_strat_risk_free)) *sqrt(252);
    sharpe_ratio = mean(pl_strat) / std(pl_strat) * sqrt(252);
    
    len = length(spx_portfolio);
    pl_strat = pl_strat';
    strat_porfolio = strat_porfolio';
    plot([spx_portfolio strat_porfolio]);
    grid on;
    
    strat_returns = Compute_Returns(strat_porfolio);
    %strat_returns(strat_returns==0) = [];
    
    %Compute yearly percentages for the strategy
    if(length(strat_porfolio) > 3000)
        tmp = strat_porfolio(1:365:len);
        percentages = diff(tmp)./tmp(1:end-1);
        percentages = percentages + 1;

        tmp = pl_spx_cum(1:365:len);
        percentages2 = diff(tmp)./tmp(1:end-1);
        percentages2 = percentages2 + 1;
    
        Y = [percentages-1, percentages2-1, repmat(risk_free_rate-1, length(percentages), 1)];
        figure;
        barh(Y);
        legend('strat', 'spx', 'risk free rate');
        
        annualized_vol = std(percentages);
        annualized_ret = mean(percentages)-1;
        
    else
        percentages = 0;
        annualized_vol = std(strat_returns);
        annualized_ret = mean(strat_returns);
    end
    

    correlation_with_market = corr(pl_strat, pl_spx);
 
    MDD_strat = maxdrawdown(strat_porfolio);
    MDD_spx = maxdrawdown(spx_portfolio);
    
    max_daily_ret = max(strat_returns);
    min_daily_ret = min(strat_returns);
    
    skew = skewness(strat_returns);
    kurt = kurtosis(strat_returns);
    
    cumulative_profit = (strat_porfolio(end) / initial_bet - 1)*100;
    
    s = struct('skew', skew, 'kurt', kurt, 'sharpe_free', sharpe_ratio_free, 'sharpe_absolute', sharpe_ratio, 'percentages', percentages, 'correlation_with_market', correlation_with_market, 'MDD_strat', MDD_strat, 'MDD_spx', MDD_spx, 'annualized_vol', annualized_vol, 'annualized_ret', annualized_ret, 'max_daily_ret', max_daily_ret, 'min_daily_ret', min_daily_ret, 'cumulative_profit', cumulative_profit, 'strat_returns', strat_returns);
    
end