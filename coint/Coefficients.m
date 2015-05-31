%Best tuple
% JPMORGAN CHASE & CO	713	 Financials	 
%PNC FINANCIAL SERVICES GROUP	949	 Financials	 
%WELLS FARGO & CO	1187	 Financials	
%p value r0 = 0.004559	
% pvalue r1 = 0.18676	
% pvalue r2 = 0.189698	
% corr R1 R2 = 0.899406	
% corr R2 R3 = 0.92561	
% corr R1 R3 = 0.836641	
% corr ret R1 R2 = 0.675507	
% corr ret R2 R3 = 0.709919	
% corr ret R1 R3 = 0.681343

clear;
load spx.mat;

jpm = Convert(pp, 713);
pnc = Convert(pp, 949);
fargo = Convert(pp, 1187);

Y = [jpm, pnc, fargo];

tbl = table(jpm, pnc, fargo);
mod = fitlm(tbl, 'jpm ~ pnc + fargo - 1');
beta = mod.Coefficients.Estimate;

%         Estimate       SE        tStat       pValue  
%pnc      0.62883     0.0075108    83.724             0
%fargo    0.14872       0.01611    9.2316    3.6069e-20

spr = jpm - beta(1) * pnc - beta(2) * fargo;
[~, pValue] = adftest(spr, 'model','TS','lags', 0);
pValue
%stationary
%[~, pValue] = adftest(cumsum(randn(1000,1)), 'model','TS','lags', 0); pValue

beg = 1000;
dur = 1000;
%spr = spr(beg:beg+dur);
[mid, uppr, lowr] = bollinger(spr);
plot([mid, uppr, lowr, spr]);

freq = Revertion_Frequency(spr, 20);
freq