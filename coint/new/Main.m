run('Init.m');
spreads = SpreadFinder(pp, 'hello.txt', 0.8);

for spread = spreads
    boll_conf = struct('wsize', 74, 'wts', 1, 'nstd', 2.5);
    [pl,balance_cum] = SimpleTradingStrategy( pp, spread, 1, 6100, boll_conf );
    plot(balance_cum);

    PerformanceAssessment(balance_cum, spx);
end