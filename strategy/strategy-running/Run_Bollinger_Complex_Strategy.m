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
% Run Simple before and load spreads_ids and spreads
%%%%%%%%%%%%%%

N = length(spreads(1).px)+100;
for spr_id = spreads_ids
    spread      = spreads(spr_id);
    [ pmcmc ]   = Sto_Volatility_Estimator( spread.px, 5000, N );
    fname       = sprintf('pmcmc_%d.mat', spr_id); 
    save(fname, 'pmcmc');
end
