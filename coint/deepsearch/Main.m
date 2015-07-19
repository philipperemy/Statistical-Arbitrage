run('../Init.m');
run('Headers.m');
spreads = SpreadFinder(pp, 0.8, 236:240); %Why there is only 7 spreads?
%range 236:240 is enough
tic;
nstd_range = 1.0:0.1:3.0;
wts_range = 0:1;
wsize_range = 5:1:100;

c = 1;
spread_length = 3;
for spread = spreads(1:spread_length)
    [ profits, vols, sharpes, I ] = CrossValidationParallel( pp, spread, @SimpleTradingStrategy, nstd_range, wts_range, wsize_range );
    I(:,4)          = sharpes;
    I(:,5)          = profits;
    I(:,6)          = vols;
    
    I               = I(all(~isnan(I),2),:); % for nan - rows
    I               = sortrows(I, 4);
    
    wsize           = I(end,3);
    wts             = I(end,2);
    nstd            = I(end,1);
    [~, balance_cum]= SimpleTradingStrategy( pp, spread, 1, length(spread.px), struct('wsize', wsize, 'wts', wts, 'nstd', nstd), 0 );
    balances(c,:)   = balance_cum;
    fprintf('spread %i analysed\n', c);
    c = c + 1;
end

w = repmat(1/spread_length, spread_length, 1);
portfolio = balances'*w;
plot(portfolio);
pl = diff(portfolio);
sharpe = (mean(pl)/std(pl))*sqrt(252);
sharpe
toc;

[sharpe_spx, sharpe_free, sharpe_absolute, percentages ] = PerformanceAssessment( portfolio, spx );