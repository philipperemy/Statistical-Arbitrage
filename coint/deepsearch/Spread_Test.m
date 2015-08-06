load Data_Canada
Y = Data(:,3:end);
[ res, spread, pvals_jci ] = SpreadConstructor([1 2 3], Y);
assert(res == 1);