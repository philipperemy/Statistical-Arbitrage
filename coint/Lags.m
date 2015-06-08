clear;
load spx.mat;
[ pValue, lm_formula, beta, spr ] = AugmentedDickeyFullerTest( pp, 236, 315, 349 );

%Can be different with regard to the lags
[h,pValue,stat,cValue,reg] = adftest(spr, 'model','TS','lags', 1:6);
h

%https://tel.archives-ouvertes.fr/hal-01020405/document
max_lag = 12*(length(Convert(pp, 236))/100)^0.25;
for p = 0:max_lag
    [h,pValue,stat,cValue,reg] = adftest(spr, 'model','TS','lags', p);
    aic(p+1) = reg.AIC;
    bic(p+1) = reg.BIC;
    ll(p+1) = reg.LL;
    hqc(p+1) = reg.HQC;
    stats(p+1) = stat;
    fprintf('computing for lag %i\n', p);
end

find(stats == min(stats))
find(bic == min(bic))
find(ll == min(ll))
find(hqc == min(hqc))
find(aic == min(aic))

subplot(3,2,1);
plot(aic);
title('aic');
subplot(3,2,2);
plot(bic);
title('bic');
subplot(3,2,3); 
plot(ll);
title('ll');
subplot(3,2,4);
plot(hqc);
title('hqc');
subplot(3,2,5);
plot(stats);
title('stats');

%%%% Other spread %%%%
[ pValue, lm_formula, beta, spr2 ] = AugmentedDickeyFullerTest( pp, 713, 949, 1187 );
max_lag = 12*(length(Convert(pp, 236))/100)^0.25;
[h,pValue,stat,cValue,reg] = adftest(spr2, 'model','TS','lags', 1:33);

for p = 0:max_lag
    [h,pValue,stat,cValue,reg] = adftest(spr2, 'model','TS','lags', p);
    aic(p+1) = reg.AIC;
    bic(p+1) = reg.BIC;
    ll(p+1) = reg.LL;
    hqc(p+1) = reg.HQC;
    stats(p+1) = stat;
    fprintf('computing for lag %i\n', p);
end

find(stats == min(stats))
find(bic == min(bic))
find(ll == min(ll))
find(hqc == min(hqc))

%The idea is that spr2 is more stable than spr
pptest(spr, 'model', 'TS', 'lags', 1:max_lag)
pptest(spr2, 'model', 'TS', 'lags', 1:max_lag)

[h,pValue,stat,cValue,mles] = jcitest([Convert(pp, 713), Convert(pp, 949), Convert(pp, 1187)], 'lags', 0:30)
[h,pValue,stat,cValue,mles] = jcitest([Convert(pp, 236), Convert(pp, 315), Convert(pp, 349)], 'lags', 0:30)

[ pValue, lm_formula, beta, spr ] = AugmentedDickeyFullerTest( pp, 245, 417, 691 );
