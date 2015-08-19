function [ s ] = PerformanceAssessment( pl_strat_cum, spx, initial_bet )

    T_strat = length(pl_strat_cum);
    T_spx = length(spx);
    T_min = min(T_spx, T_strat);
    pl_strat_cum = pl_strat_cum(1:T_min);
    spx = spx(1:T_min);

    if(~exist('initial_bet', 'var'))
        initial_bet = 10000;
    end
    
    percent_diff_spx = diff(spx)./spx(1:end-1,:);
    pl_spx_cum = initial_bet*cumprod([0 percent_diff_spx']+1)';
    pl_spx = diff(pl_spx_cum);

    diff_pl_strat_cum = diff(pl_strat_cum);

    %if the risk free rate is 2% = *1.02 per year
    risk_free_rate = 1.02;
    daily_percentage = risk_free_rate^(1/365);
    daily_risk_free_rate = repmat(daily_percentage, 1, length(pl_strat_cum));
    pl_risk_free_cum = initial_bet*cumprod(daily_risk_free_rate)';
    pl_risk_free = diff(pl_risk_free_cum);
   
    norm_pl_rf = diff_pl_strat_cum - mean(pl_risk_free);
    sharpe_free = (mean(norm_pl_rf) / std(norm_pl_rf))*sqrt(252);

    sharpe_absolute = mean(diff_pl_strat_cum) / std(diff_pl_strat_cum) * sqrt(252);
    
    len = length(pl_spx_cum);
    flat_init_bet_cum = repmat(initial_bet, len, 1);

    diff_pl_strat_cum = diff_pl_strat_cum';
    pl_strat_cum = pl_strat_cum';
    plot([pl_risk_free_cum pl_spx_cum pl_strat_cum flat_init_bet_cum]);
    
    strat_returns = Compute_Returns(pl_strat_cum);
    strat_returns(strat_returns==0) = [];
    
    %Compute yearly percentages for the strategy
    if(T_min > 3000)
        tmp = pl_strat_cum(1:365:len);
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
        percentages2 = 0;
        annualized_vol = std(strat_returns);
        annualized_ret = mean(strat_returns);
    end
    
    correlation_with_market = corr(diff_pl_strat_cum, pl_spx);
    
    
    MDD_strat = maxdrawdown(pl_strat_cum);
    MDD_spx = maxdrawdown(spx);
    
    max_daily_ret = max(strat_returns);
    min_daily_ret = min(strat_returns);
    
    skew = skewness(strat_returns);
    kurt = kurtosis(strat_returns);
    
    cumulative_profit = (pl_strat_cum(end) / initial_bet - 1)*100;
    
    s = struct('skew', skew, 'kurt', kurt, 'sharpe_free', sharpe_free, 'sharpe_absolute', sharpe_absolute, 'percentages', percentages, 'correlation_with_market', correlation_with_market, 'MDD_strat', MDD_strat, 'MDD_spx', MDD_spx, 'annualized_vol', annualized_vol, 'annualized_ret', annualized_ret, 'max_daily_ret', max_daily_ret, 'min_daily_ret', min_daily_ret, 'cumulative_profit', cumulative_profit, 'strat_returns', strat_returns);
    
end


%     trading costs? Costs cannot be handled here. Or it's more difficult.
%     costs = 1.01;
%     trades = sum(pl ~= 0);
%     %pl_strat_cum = cumsum(pl)+initial_bet;
%     pl_strat_cum = cumsum(pl*(initial_bet/100))+initial_bet;
%     trades_outstanding = (pl_strat_cum & (pl~=0)).*pl_strat_cum;
%     %costs_outstanding = -trades_outstanding*(costs-1); %Here you trade always with ONE share. You don't invest money again.
%     costs_outstanding = cumsum(-(pl ~= 0));
%     plot([(pl_strat_cum + costs_outstanding)' pl_strat_cum']);

%     pl_strat_w_cost = diff(pl_strat_cum + costs_outstanding);
%     norm_pl_w_cost_rf = pl_strat_w_cost - mean(pl_risk_free);
%     sharpe_free_w_cost = (mean(norm_pl_w_cost_rf) / std(norm_pl_w_cost_rf))*sqrt(252);
%     sharpe_free_w_cost

%     b = pl_spx - mean(pl_risk_free);
%     sharpe_spx_risk_free = mean(b)/std(b)*sqrt(252);
%     sharpe_spx_risk_free
