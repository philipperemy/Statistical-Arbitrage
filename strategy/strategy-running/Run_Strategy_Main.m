clear;clc;

addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('../../pmcmc');
addpath('../../filters');
addpath('../../models');
addpath('..');
load ../../data/spx.mat;
format longg;


ZSCORE_ID           = 3;
BOLLINGER_CPX       = 2;
BOLLINGER_SIMPLE    = 1;

strat_id            = BOLLINGER_CPX;

% BOLLINGER_SIMPLE = 0.4614 
% BOLLINGER_CPX = 0.4321,
% load('spreads_1_731.mat');
% load('spx_1_731.mat');
% [ portfolio_cumsum_1_731 ] = Run_Strategy( pp,spreads, spx, BOLLINGER_SIMPLE );
% fprintf('done 1\n');
% 
% load('spreads_731_1462.mat');
% load('spx_731_1462.mat');
% [ portfolio_cumsum_731_1462 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 731\n');
% 
% load('spreads_1462_2192.mat');
% load('spx_1462_2192.mat');
% [ portfolio_cumsum_1462_2192 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 1462\n');
% 
% load('spreads_2192_2923.mat');
% load('spx_2192_2923.mat');
% [ portfolio_cumsum_2192_2923 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 2192\n');
% 
% load('spreads_2923_3653.mat');
% load('spx_2923_3653.mat');
% [ portfolio_cumsum_2923_3653 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 2923\n');
% 
% load('spreads_3653_4384.mat');
% load('spx_3653_4384.mat');
% [ portfolio_cumsum_3653_4384 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 3653\n');
% 
% load('spreads_4384_5114.mat');
% load('spx_4384_5114.mat');
% [ portfolio_cumsum_4384_5114 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 4384\n');
% 
% load('spreads_5114_5845.mat');
% load('spx_5114_5845.mat');
% [ portfolio_cumsum_5114_5845 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 5114\n');
% 
% load('spreads_5845_6575.mat');
% load('spx_5845_6575.mat');
% [ portfolio_cumsum_5845_6575 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 5845\n');

% Line 438: i = 20, CR = 1.277284e+04
% Line 444: i = 15, CR = 1.087485e+04
% Line 448: i = 34, CR = 1.011817e+04
% Line 485: i = 28, CR = 2.873406e+03
% Line 555: i = 24, CR = 1.253293e+05
% Line 603: i = 10, CR = 9.018609e+03
% Line 636: i = 31, CR = 9.112141e+03
% Line 700: i = 5, CR = 6.868178e+04
% Line 840: i = 14, CR = 1.596698e+04
% Line 893: i = 27, CR = 1.026277e+05
% Line 940: i = 19, CR = 1.289233e+04
% Line 1065: i = 4, CR = 1.018068e+04
% Line 1074: i = 30, CR = 1.211030e+04
% Line 1114: i = 23, CR = 1.253293e+05
% Line 1116: i = 9, CR = 8.768622e+03
% Line 1273: i = 26, CR = 9.822859e+03
% Line 1408: i = 13, CR = 1.578188e+04
% Line 1433: i = 18, CR = 2.550341e+04
% Line 1442: i = 29, CR = 8.419956e+03
% Line 1530: i = 8, CR = 1.056858e+04
% Line 1668: i = 22, CR = 9.651261e+03
% Line 1766: i = 25, CR = 1.152300e+04
% Line 1801: i = 12, CR = 8.629521e+03
% Line 1866: i = 3, CR = 2.120840e+04
% Line 1890: i = 17, CR = 1.724399e+04
% Line 1911: i = 36, CR = 5.102085e+03
% Line 1924: i = 7, CR = 1.183115e+04
% Line 2122: i = 21, CR = 9.515091e+03
% Line 2158: i = 38, CR = 1.666710e+04
% Line 2271: i = 11, CR = 1.487422e+04
% Line 2273: i = 35, CR = 1.170330e+04
% Line 2373: i = 2, CR = 1.782267e+04
% Line 2401: i = 33, CR = 7.145196e+03
% Line 2427: i = 16, CR = 1.069201e+04
% Line 2442: i = 6, CR = 1.498092e+04
% Line 2672: i = 40, CR = 1.050681e+04
% Line 2676: i = 42, CR = 6.160043e+03
% Line 2682: i = 43, CR = 1.091262e+04
% Line 2800: i = 1, CR = 1.465440e+04
% Line 2832: i = 32, CR = 1.185243e+04
% Line 2877: i = 44, CR = 2.304224e+04
% Line 2955: i = 45, CR = 1.919795e+04
% Line 3078: i = 39, CR = 1.219933e+04
% Line 3139: i = 41, CR = 1.042181e+04
% Line 3154: i = 47, CR = 1.872888e+04
% Line 3164: i = 49, CR = 1.950518e+04
% Line 3192: i = 48, CR = 1.661736e+04
% Line 3283: i = 46, CR = 1.228488e+04
% Line 3288: i = 37, CR = 1.130668e+04

% mean(sr)
% 
% ans =
%          0.419872057980445
% WITH SIMPLE IT WAS 0.369761053758438

load('spreads_6575_7306.mat');
load('spx_6575_7306.mat');
[ portfolio_cumsum_6575_7306 ] = Run_Strategy( pp,spreads, spx, 1 );
fprintf('done 6575\n');
% 
% load('spreads_7306_8036.mat');
% load('spx_7306_8036.mat');
% [ portfolio_cumsum_7306_8036 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 7306\n');
% 
% load('spreads_8036_8767.mat');
% load('spx_8036_8767.mat');
% [ portfolio_cumsum_8036_8767 ] = Run_Strategy( pp,spreads, spx, strat_id );
% fprintf('done 8036\n');