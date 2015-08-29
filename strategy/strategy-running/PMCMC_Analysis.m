addpath('../../helpers/');
addpath('../../coint/deepsearch');
addpath('../../coint/impl');
addpath('../../pmcmc');
addpath('../../filters');
addpath('../../models');
addpath('..');
load ../../data/spx.mat;
format longg;

%%% LOAD SPREAD MODULE %%%
load('spreads_4384_5114.mat');

initial_bet  = 10000;
spread_count = length(spreads);
mat_tr 		 = zeros(spread_count , 5 );
boll_conf    = struct('wsize', 20, 'wts', 1, 'nstd', 2);
zscore_conf  = struct('sell_open', 1, 'sell_close', -1, 'buy_open', -1, 'buy_close', 1);

for i = 1:spread_count
    try
        fprintf('TRAIN Bollinger Complex, i = %d\n', i);
        s = sprintf('pmcmc_%d.mat', i);
        load(s);
        mat = complex.mat;
        complex_st(i) = complex;
        mat_tr(i,1) = mat(1);
        mat_tr(i,2) = mat(2);
        mat_tr(i,3) = mat(3);
        mat_tr(i,4) = mat(4);
        mat_tr(i,5) = mat(5);
    catch
    end
end

N          = length(complex_st);
theta_mat  = zeros(1,6);
for i = 1:N
   cpx_st = complex_st(i);
   try
       rho1 = cpx_st.model.pmcmc.rho1;
       sigma1 = cpx_st.model.pmcmc.sigma1;
       rho2 = cpx_st.model.pmcmc.rho2;
       sigma2 = cpx_st.model.pmcmc.sigma2;
       beta = cpx_st.model.pmcmc.beta;
       cor = cpx_st.model.pmcmc.cor;
       theta_mat(end+1,:) = [rho1 sigma1 rho2 sigma2 beta cor];
   catch
   end
end

rho_props = theta_mat(:,1);
rho_props(rho_props == 0) = [];