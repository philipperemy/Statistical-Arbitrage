clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');
addpath('./surf/');
addpath('../coint/deepsearch/');
load('../data/spx.mat');

pp = TruncateData(pp, 5000, 7000);
load('spreads.mat');
nstd_range = 1.0:0.1:3.0;
wsize_range = 5:1:100;
wts = 0;

% [ spreads ] = SpreadFinder( pp, 0.8, 460:480);
% SpreadProcessor(spreads) 
% SpreadProcessor ( spreads );

step = 0;
Z_cumsum = zeros(96,21); %same range
for i = 1:50:length(spreads)
    Spread = spreads(i);
    [ profits, vols, sharpes, I ] = CrossValidationParallel( pp, Spread, @SimpleTradingStrategy, nstd_range, wts, wsize_range );
    I(:,4)          = sharpes;
    I(:,5)          = profits;
    I(:,6)          = vols;

    I(isnan(I)) = 0;
    A = [I(:,3) I(:,1) I(:,4)];
    W = sortrows(A, 3);

    wsize           = W(end,1);
    nstd            = W(end,2);
    [~, balance_cum]= SimpleTradingStrategy( pp, Spread, 1, length(Spread.px), struct('wsize', wsize, 'wts', wts, 'nstd', nstd), 0 );

    close all;
    figure;
    Z = Draw_Surf(I, nstd_range, wsize_range);
    drawnow;
    fprintf('id = %i, balance = %f, sharpe = %f\n', i, balance_cum(end), W(end,3));
    Z_cumsum = Z_cumsum + Z;
    step = step+1;
end
