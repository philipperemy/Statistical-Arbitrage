clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../models/');

load('../data/spx.mat');

p_y_given_x = @(y,sig) normpdf(y, 0, sig);

APPLE_TICKER = 239;
appl_prices = GetPrice(pp, APPLE_TICKER);

log_returns = diff(log(appl_prices))';
log_returns = log_returns(4000:5000);

beta = 0.5;
st = SVLogReturns(log_returns, beta);

p_y_given_th = BootstrapParticleFilter(st.y, 0.99, 0.15, 0.5, 100000, p_y_given_x);
true_val = 2.1694e+03;
conf_int = [0.99*true_val 1.01*true_val];
assert(conf_int(1) < p_y_given_th & p_y_given_th < conf_int(2));

p_y_given_th = BootstrapParticleFilter(st.y, 0.99, 0.15, 0.5, 10000, p_y_given_x);
true_val = 2.1694e+03;
conf_int = [0.99*true_val 1.01*true_val];
assert(conf_int(1) < p_y_given_th & p_y_given_th < conf_int(2));

p_y_beta_1 = BootstrapParticleFilter(st.y, 0.99, 0.15, 0.01, 10000, p_y_given_x);
p_y_beta_2 = BootstrapParticleFilter(st.y, 0.99, 0.15, 0.5, 10000, p_y_given_x);
assert(p_y_beta_2 < p_y_beta_1);

%test stability.
N_seq = 100:100:5000;
MC_steps = 50;
p_y_given_th_mat = zeros(length(N_seq), MC_steps);
i = 1;
for N = N_seq
    fprintf('N = %i\n', N);
    parfor M = 1:MC_steps
        %mean to run parallel and to smooth the results
        %can be parallelized easily
        a1 = BootstrapParticleFilter(st.y, 0.99, 0.15, st.beta, N, p_y_given_x);
        a2 = BootstrapParticleFilter(st.y, 0.99, 0.15, st.beta, N, p_y_given_x);
        a3 = BootstrapParticleFilter(st.y, 0.99, 0.15, st.beta, N, p_y_given_x);
        a4 = BootstrapParticleFilter(st.y, 0.99, 0.15, st.beta, N, p_y_given_x);
       p_y_given_th_mat(i,M) = mean([a1 a2 a3 a4]);
    end
    i = i + 1;
end
vars = var(p_y_given_th_mat, 0,2);
assert(any(vars(20:end)' < 10));

%%%%%
% Idea bollinger bands for SV model
%%%%%

% est_X = zeros(10, length(st.y));
% for i = 1:10
%     [p_y_given_theta, estimated_states] = BootstrapParticleFilter(st.y, 0.99, 0.15, 0.5, 10000, p_y_given_x);
%     est_X(i,:) = estimated_states;
% end
% vt = mean(est_X,1);
% exp(normrnd(0, 0.5*exp(vt/2)));
% 
% px = appl_prices(4000:5000)';
% px_lag = lagmatrix(px, 1);
% 
% px_mat = zeros(10, length(px_lag));
% for i = 1:10
%    px_mat(i,:) = px_lag.*exp(normrnd(0, 0.5*exp(vt/2)))'; 
% end
% plot([quantile(px_mat, 0.9, 1)' px' quantile(px_mat, 0.1, 1)']);
%E[x_t|x_t-1] because x_t is based on y_t which is not measurable wrt Ft

fprintf('Test completed\n');