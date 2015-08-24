clear;clc;

ZSCORE_ID           = 3;
BOLLINGER_CPX       = 2;
BOLLINGER_SIMPLE    = 1;

strat_id            = BOLLINGER_CPX;

load('spreads_1_731.mat');
load('spx_1_731.mat');
[ portfolio_cumsum_1_731 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 1\n');

load('spreads_731_1462.mat');
load('spx_731_1462.mat');
[ portfolio_cumsum_731_1462 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 731\n');

load('spreads_1462_2192.mat');
load('spx_1462_2192.mat');
[ portfolio_cumsum_1462_2192 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 1462\n');

load('spreads_2192_2923.mat');
load('spx_2192_2923.mat');
[ portfolio_cumsum_2192_2923 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 2192\n');

load('spreads_2923_3653.mat');
load('spx_2923_3653.mat');
[ portfolio_cumsum_2923_3653 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 2923\n');

load('spreads_3653_4384.mat');
load('spx_3653_4384.mat');
[ portfolio_cumsum_3653_4384 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 3653\n');

load('spreads_4384_5114.mat');
load('spx_4384_5114.mat');
[ portfolio_cumsum_4384_5114 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 4384\n');

load('spreads_5114_5845.mat');
load('spx_5114_5845.mat');
[ portfolio_cumsum_5114_5845 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 5114\n');

load('spreads_5845_6575.mat');
load('spx_5845_6575.mat');
[ portfolio_cumsum_5845_6575 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 5845\n');

load('spreads_6575_7306.mat');
load('spx_6575_7306.mat');
[ portfolio_cumsum_6575_7306 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 6575\n');

load('spreads_7306_8036.mat');
load('spx_7306_8036.mat');
[ portfolio_cumsum_7306_8036 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 7306\n');

load('spreads_8036_8767.mat');
load('spx_8036_8767.mat');
[ portfolio_cumsum_8036_8767 ] = Run_Strategy( spreads, spx, strat_id );
fprintf('done 8036\n');