function [sharpe_spx, sharpe_free, sharpe_absolute, percentages ] = PerformanceAssessment( pl_strat_cum, spx, initial_bet )

    if(~exist('initial_bet', 'var'))
        initial_bet = 10000;
    end
    
    percent_diff_spx = diff(spx)./spx(1:end-1,:);
    pl_spx_cum = initial_bet*cumprod(percent_diff_spx+1);
    pl_spx = diff(pl_spx_cum);

    norm_pl = diff(pl_strat_cum) - mean(pl_spx); %and SPX is not a risk free!
    sharpe_spx = (mean(norm_pl) / std(norm_pl))*sqrt(252);

    %if the risk free rate is 2% = *1.02 per year
    risk_free_rate = 1.02;
    daily_percentage = risk_free_rate^(1/365);
    daily_risk_free_rate = repmat(daily_percentage, 1, length(pl_strat_cum));
    pl_risk_free_cum = initial_bet*cumprod(daily_risk_free_rate)';
    pl_risk_free = diff(pl_risk_free_cum);

    norm_pl_rf = diff(pl_strat_cum) - mean(pl_risk_free);
    sharpe_free = (mean(norm_pl_rf) / std(norm_pl_rf))*sqrt(252);

    sharpe_absolute = mean(diff(pl_strat_cum)) / std(diff(pl_strat_cum)) * sqrt(252);
    
    len = length(pl_spx_cum);
    flat_init_bet_cum = repmat(initial_bet, len, 1);

    plot([pl_risk_free_cum pl_spx_cum pl_strat_cum flat_init_bet_cum]);
    figure;
    
    %Compute yearly percentages for the strategy
    tmp = pl_strat_cum(1:365:len);
    percentages = diff(tmp)./tmp(1:end-1);
    percentages = percentages + 1;

    tmp = pl_spx_cum(1:365:len);
    percentages2 = diff(tmp)./tmp(1:end-1);
    percentages2 = percentages2 + 1;

    Y = [percentages-1, percentages2-1, repmat(risk_free_rate-1, length(percentages), 1)];
    barh(Y);
    legend('strat', 'spx', 'risk free rate');
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
