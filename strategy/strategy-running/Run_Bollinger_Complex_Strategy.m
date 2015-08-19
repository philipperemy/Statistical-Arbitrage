addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('../../pmcmc');
addpath('../../filters');
addpath('../../models');
addpath('..');
load ../../data/spx.mat;
format longg;

%%%%%%%%%%%%%%
% Run Simple before and load Spreads_PMCMC
%%%%%%%%%%%%%%

for i = 1:length(Spreads_PMCMC)
    spread = Spreads_PMCMC(i);
    [ vol_two_factors_leverage ] = Sto_Volatility_Estimator( spread.px, 501, 1000 );
end
