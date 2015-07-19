clear;
addpath('../filters/');
addpath('../helpers/');
addpath('../likelihoods/');
addpath('../models/');
addpath('../pmcmc/');

mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

rho = 0.91;       %Scale parameter for X
sigma = 1;        %Standard deviation of the X (process)
beta = 0.5;       %Scale parameter for Y
steps = 100;      %Number of steps
nu = 3;           %Degrees of freedom of the innovation of the measurement
st = StochasticVolatilityModelStudent(rho, sigma, beta, steps, nu);

N = 200;
MAX = 50;
rho_seq = -0.99:0.01:0.99;
log_marginal_likelihood_mat = zeros(length(rho_seq), MAX);
j = 1;
for rho = rho_seq
    for i = 1:MAX
        log_marginal_likelihood_mat(j, i) = Ex4_Compute_Log_Likelihood(st.y, rho, st.sigma, st.beta, N, st.nu);
        fprintf('%i, %i\n', j, i);
    end
    j = j + 1;  
end

plot(log_marginal_likelihood_mat); %Constant, at least at the equilibrium

%sigma = 3, nu = 6.
j = 1;
log_marginal_likelihood_mat2 = zeros(length(rho_seq), MAX);
for rho = rho_seq
    for i = 1:MAX
        log_marginal_likelihood_mat2(j, i) = Ex4_Compute_Log_Likelihood(st.y, rho, 3, st.beta, N, 6);
        fprintf('%i, %i\n', j, i);
    end
    j = j + 1;  
end

plot(std(log_marginal_likelihood_mat2, 0, 2));

%Changing the number of particles (equilibrium)
N_seq = 10:10:1000;
j = 1;
log_marginal_likelihood_mat3 = zeros(length(N_seq), MAX);
for N_it = N_seq
    for i = 1:MAX
        log_marginal_likelihood_mat3(j, i) = Ex4_Compute_Log_Likelihood(st.y, st.rho, st.sigma, st.beta, N_it, st.nu);
        fprintf('%i, %i\n', j, i);
    end
    j = j + 1;
end
plot(log_marginal_likelihood_mat3);

j = 1;
log_marginal_likelihood_mat4 = zeros(length(N_seq), MAX);
for N_it = N_seq
    for i = 1:MAX
        log_marginal_likelihood_mat4(j, i) = Ex4_Compute_Log_Likelihood(st.y, st.rho, 4, st.beta, N_it, 12);
        fprintf('%i, %i\n', j, i);
    end
    j = j + 1;
end
plot(log_marginal_likelihood_mat4);

%When std goes to zero, it means the convergence is done.
STD = [std(log_marginal_likelihood_mat3,0,2) std(log_marginal_likelihood_mat4,0,2)];
plot(N_seq, STD);

%unbiased
a = mean(log_marginal_likelihood_mat3, 2);
true_val = a(end);
clear a;

plot(N_seq, mean(abs(log_marginal_likelihood_mat3 - true_val),2));
bar(N_seq(2:end), diff(mean(log_marginal_likelihood_mat3,2)));

quantile(STD(:,1), 0.05)
quantile(STD(:,2), 0.05)

ideal_particles = find(STD < 1 == 1);
N_ideal = ideal_particles(1); %%max or mean.
N_ideal = N_seq(N_ideal);
clear ideal_particles;

%Now bruteforce it.
%Changing the number of particles (equilibrium)
%CAN BE USED FOR A FIGURE! VERY GOOD. BELOW 1. WHAT WE SHOULD HAVE

%%%%%%%%%%%%%% THE INTERESTING PART STARTS HERE %%%%%%%%%%%%%%
N_seq = 100:100:3000;
MAX = 100;
clc;
rho = 0.91;       %Scale parameter for X
sigma = 1;        %Standard deviation of the X (process)
beta = 0.5;       %Scale parameter for Y
steps = 1000;      %Number of steps
nu = 3;           %Degrees of freedom of the innovation of the measurement
st = StochasticVolatilityModelStudent(rho, sigma, beta, steps, nu);
particles_dist_test = ParticleCountProfilingStudentSV_Unit( st, N_seq, 100, rho, sigma, beta, nu );
hold off
plot(100:100:3000, var(particles_dist_test, 0, 2))
hold on
plot([1 3000],[0.8464 0.8464])
legend('Variance of p_N(y|\theta)', 'Optimal variance');


log_marginal_likelihood_mat2 = Sensibility_with_rho( 1000, st, sigma, beta, nu );
plot(rho_seq, log_marginal_likelihood_mat2); %beautiful

plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat2, 0, 2),'s',20,1)) %variance quite stable
legend('Smoothed Variance of p_N(y|\theta)');

log_marginal_likelihood_mat3 = Sensibility_with_rho( 2000, st, sigma, beta, nu );
hold on;
plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat3, 0, 2),'s',20,1)) %variance quite stable

log_marginal_likelihood_mat4 = Sensibility_with_rho( 500, st, sigma, beta, nu );

hold on;
plot(rho_seq, tsmovavg(var(log_marginal_likelihood_mat4, 0, 2),'s',20,1)) %variance quite stable

legend('N = 1000', 'N = 2000', 'N = 500');
