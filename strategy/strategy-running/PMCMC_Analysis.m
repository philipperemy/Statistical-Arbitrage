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
rho_props(find(diff(rho_props) == 0)) = [];

sigma_props = theta_mat(:,2);
sigma_props(sigma_props == 0) = [];
sigma_props(find(diff(sigma_props) == 0)) = [];

rho2_props = theta_mat(:,3);
rho2_props(rho2_props == 0) = [];
rho2_props(find(diff(rho2_props) == 0)) = [];

sigma2_props = theta_mat(:,4);
sigma2_props(sigma2_props == 0) = [];
sigma2_props(find(diff(sigma2_props) == 0)) = [];

beta_props = theta_mat(:,5);
beta_props(beta_props == 0) = [];
beta_props(find(diff(beta_props) == 0)) = [];

cor_props = theta_mat(:,6);
cor_props(cor_props == 0) = [];
cor_props(find(diff(cor_props) == 0)) = [];

close all;
Draw_Histogram(rho_props, [0.90 1], 60, '\phi_X', 'north');
Draw_Histogram(sigma_props, [0 2], 60, '\sigma_X', 'north');
Draw_Histogram(rho2_props, [0 1], 60, '\phi_Z', 'north');
Draw_Histogram(sigma2_props, [0 2], 60, '\sigma_Z', 'north');
Draw_Histogram(beta_props, [0 4.8], 60, '\beta', 'north');
Draw_Histogram(cor_props, [-1 1], 60, '\rho', 'north');

figure;
%For the statistics
MCMC_Checks(rho_props);
MCMC_Checks(sigma_props);
MCMC_Checks(rho2_props);
MCMC_Checks(sigma2_props);
MCMC_Checks(beta_props);
MCMC_Checks(cor_props);

sr = mat_tr(:,3);
sr(sr == 0) = [];
sr(isnan(sr)) = [];
mean(sr)

Draw_Histogram(sr, [-5 5], 20, 'SR (out-sample)', 'northwest');
figure;
MCMC_Checks(sr);

%sr(find(diff(sr) == 0)) = [];
%CAN IMPROVE