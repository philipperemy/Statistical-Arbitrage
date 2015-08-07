load Data_Canada
Y = Data(:,3:end);
[ res, spread, pvals_jci ] = SpreadConstructor([1 2 3], Y);
assert(res == 1);

[h,pValue,stat,cValue,mles] = jcitest(Y,'model','H1');
YLag = Y(2:end,:);
T = size(YLag,1);
B = mles.r2.paramVals.B;
c0 = mles.r2.paramVals.c0;
plot(dates(2:end),YLag*B+repmat(c0',T,1))
grid on