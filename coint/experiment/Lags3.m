clear;
load spx.mat;
[ pValue, lm_formula, beta, spr ] = AugmentedDickeyFullerTest( pp, 236, 315, 349 );
%https://tel.archives-ouvertes.fr/hal-01020405/document
max_lag = 12*(length(Convert(pp, 236))/100)^0.25;

for p = 0:max_lag
    [h,pValue,stat,cValue,reg] = pptest(spr, 'model', 'TS', 'lags', p);
    aic_v(p+1) = reg.AIC;
    bic_v(p+1) = reg.BIC;
    ll_v(p+1) = reg.LL;
    hqc_v(p+1) = reg.HQC;
    stats_v(p+1) = stat;
    fprintf('computing for lag %i\n', p);
end

% PP TEST doc very interesting
%     The value of 'model' is determined by the growth characteristics of
%     the time series being tested, and should be chosen with a specific
%     testing strategy in mind. As discussed in Elder & Kennedy [2],
%     including too many regressors results in lost power, while including
%     too few biases the test in favor of the null. In general, if a series
%     is growing, the 'TS' model provides a reasonable trend-stationary
%     alternative to a unit-root process with drift. If a series is not
%     growing, 'AR' and 'ARD' models provide reasonable stationary
%     alternatives to a unit-root process without drift. The 'ARD'
%     alternative has mean c/(1-a); the 'AR' alternative has mean 0.
