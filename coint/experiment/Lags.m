A %non stationary
Y %stationary
%http://fr.mathworks.com/help/econ/unit-root-nonstationarity.html Schwert

series = A;
lag_max = round(12*(length(series)/100)^0.25);


%A lower AIC means a model is considered to be closer to the truth
aic = zeros(1);
bic = zeros(1);
hqc = zeros(1);
for lag = fliplr(1:lag_max)
    [h,pValue,stat,cValue,reg] = adftest(series, 'model','TS','lags', lag);
    reg.tStats.pVal
    sig_coeffs = sum(reg.tStats.pVal < 0.05);
    fprintf('lag = %i, AIC = %f, BIC = %f, HQC = %f, sig = %i\n', lag, reg.AIC, ...
        reg.BIC, reg.HQC, sig_coeffs);
    aic(lag) = reg.AIC;
    bic(lag) = reg.BIC;
    hqc(lag) = reg.HQC;
end
plot(aic);

