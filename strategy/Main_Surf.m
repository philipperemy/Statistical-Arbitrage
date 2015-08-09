clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');
addpath('./surf/');

%http://www.mathworks.com/examples/econometrics/1501-test-stock-data-for-a-random-walk

load('../data/spx.mat');
stock_ids = [713 949 1187];
stocks_prices = [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ];
[ h, Spread, ~ ] = SpreadConstructor(stock_ids, [GetPrice( pp, stock_ids(1) ), GetPrice( pp, stock_ids(2) ), GetPrice( pp, stock_ids(3)) ]);

tic;
nstd_range = 1.0:0.1:3.0;
wsize_range = 5:1:100;

wts = 0;
[ profits, vols, sharpes, I ] = CrossValidationParallel( pp, Spread, @SimpleTradingStrategy, nstd_range, wts, wsize_range );
I(:,4)          = sharpes;
I(:,5)          = profits;
I(:,6)          = vols;

I(isnan(I)) = 0;
A = [I(:,3) I(:,1) I(:,4)];
W = sortrows(A, 3);

wsize           = W(end,1);
nstd            = W(end,2);
[~, balance_cum]= SimpleTradingStrategy( pp, spread, 1, length(spread.px), struct('wsize', wsize, 'wts', wts, 'nstd', nstd), 0 );

Draw_Surf(I, nstd_range, wsize_range);
