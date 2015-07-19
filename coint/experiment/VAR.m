%ALCOA INC, 236, Materials, AVERY DENNISON CORP, 315, Materials, BANK OF NEW YORK MELLON CORP, 349, Financials, 0.011761, 0.053341, 0.058756, 0.844806, 0.881141, 0.810372, 0.435591, 0.384840, 0.411050

%http://www.memoireonline.com/11/12/6430/m_Financial-development-and-economic-growth-evidence-from-Niger8.html
%http://stats.stackexchange.com/questions/21539/what-is-the-correct-procedure-to-choose-the-lag-when-performing-johansen-cointeg

Y = [Convert(pp, 236) Convert(pp, 315) Convert(pp, 349)];
p = 20;
n = size(Y,2);

%Just calibrate a VAR model here
Spec_VARMAX = vgxset('n',n,'nAR',p,'Constant',true);
[EstSpec,EstStdErrors,LLF,W] = vgxvarx(Spec_VARMAX, Y);

Spec_VARMAX.nAR = 200;
[EstSpec2,EstStdErrors2,LLF2,W2] = vgxvarx(Spec_VARMAX, Y);

%The interesting part starts here.
max_lag = 100;
for p = 1:max_lag
    Spec_VARMAX.nAR = p;
    [~,~,llf,~] = vgxvarx(Spec_VARMAX, Y);
    [aic, bic] = aicbic(llf, p+3, 6100);
    fprintf('lag = %i, aic = %f, bic = %f\n', p, aic, bic);
    aic_v(p) = aic;
    bic_v(p) = bic;
end

%h = 0 - no autocorrelation
%h = 1 - autocorrelation
%http://stats.stackexchange.com/questions/21539/what-is-the-correct-procedure-to-choose-the-lag-when-performing-johansen-cointeg
lbqtest(W(:,1)) %portmanteau Ljung-Box Q-test for residual autocorrelation
%assesses the null hypothesis that a series of residuals exhibits no autocorrelation for a fixed number of lags L, 
%against the alternative that some autocorrelation coefficient ρ(k), k = 1, ..., L, is nonzero
lbqtest(W(:,2))
lbqtest(W(:,3))

%Same but better written
max_lag = 60;
lag_LLF = zeros(1,max_lag);
for p = 0:max_lag
    Spec_VARMAX = vgxset('n',3,'nAR',p, 'Constant',true);
    [EstSpec,EstStdErrors,LLF,W] = vgxvarx(Spec_VARMAX, Y);
    lag_LLF(p+1) = LLF;
    fprintf('computing for lag %i\n', p);
end
numObs = length(Convert(pp, 236));
[aic, bic] = aicbic(lag_LLF, 1:max_lag+1, numObs);


