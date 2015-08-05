run('../Init.m');
addpath('../../helpers');

ids = [281 289 726];
[ ~, spread, ~ ] = SpreadConstructorLog(ids, log(GetPriceArray(pp, ids)));
ZScore(spread.px);