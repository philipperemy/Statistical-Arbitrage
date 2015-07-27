clear;
clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../likelihoods/');
addpath('../models/');
addpath('../pmcmc/');

mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

% [var_part_dist, N_seq] = ParticleDistribVarTest(200);
% hold off
% plot(N_seq, var_part_dist)
% hold on
% plot([1 3000],[0.8464 0.8464])
% legend('Variance of p_N(y|\theta)', 'Optimal variance');
% xlabel('Number of particles (N)');
% ylabel('Var(p_N(y|\theta))');

N_seq = 100:100:3000;
T = 100:100:2000;
len_N_seq = length(N_seq);
len_T_seq = length(T);
var_part_dist_mat = zeros(len_T_seq, len_N_seq);
c = 1;
for t = T
    fprintf('t = %i\n', t);
    [var_part_dist, N_seq] = ParticleDistribVarTest(t);
    var_part_dist_mat(c, :) = var_part_dist;
    c = c + 1;
end

X = zeros(1,1);
Y = zeros(1,1);
Z = zeros(1,1);
for c = 1:len_T_seq
    try
        id = sum(var_part_dist_mat(c, :) > 0.8484); %0.8484
        fprintf('T=%i, opt N = %i\n', T(c), N_seq(id));
        X(end+1) = T(c);
        Y(end+1) = N_seq(id);
        Z(end+1) = 1.20 * T(c) + 3; %polyfit(X,Y,1)
    catch ex
    end
end
plot(X,Y,'b--o', X, Z, 'r');
ylabel('Optimal number of particles N');
xlabel('Number of steps T');
legend('N_{opt} with var 0.8484', '1.20 \times T + 3.00', 'Location','southeast');

[log_marginal_likelihood_mat2] = Sensibility_with_rho( 1000, 0.70:0.01:0.99, 0.91, 1, 0.2, 3, 1000 );
plot(rho_seq, log_marginal_likelihood_mat2); %beautiful

plot(0.70:0.01:0.99, log_marginal_likelihood_mat2, 'b.');
ylabel('Log likelihood');
xlabel('\rho');
legend('log p_N(y|\theta)');


plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat2, 0, 2),'s',20,1)) %variance quite stable
legend('Smoothed Variance of p_N(y|\theta)');

log_marginal_likelihood_mat3 = Sensibility_with_rho( 2000, 0.91, 1, 0.2, 3, 1000 );
hold on;
plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat3, 0, 2),'s',20,1)) %variance quite stable

log_marginal_likelihood_mat4 = Sensibility_with_rho( 500, 0.91, 1, 0.2, 3, 1000 );

hold on;
plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat4, 0, 2),'s',20,1)) %variance quite stable

legend('N = 1000', 'N = 2000', 'N = 500');
