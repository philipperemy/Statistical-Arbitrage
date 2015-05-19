mySeed = 57; % an integer number
rng(mySeed,'twister') %You can replace 'twister' with other generators

clear;
rho = 0.91;       %Scale parameter for X
sigma = 1;        %Standard deviation of the X (process)
beta = 0.5;       %Scale parameter for Y
N = 1000;         %Number of particules
steps = 100;      %Number of steps
nu = 3;

x = zeros(1, steps+1);
y = zeros(1, steps);
x(1) = randn * sqrt(sigma^2/(1-rho^2));      

for i=1:steps
    x(i+1) = rho  * x(i)  + sigma*randn;
    y(i)   = beta * exp(0.5*x(i))*trnd(nu); %randn to test and nu is high
end
x = x(1:steps);

%rho and sigma are unknown
%beta is known
steps_mcmc = 100;

rho_prop = zeros(1, steps_mcmc);
sigma_prop = zeros(1, steps_mcmc);
nu_prop = zeros(1, steps_mcmc);

%http://www.mathworks.com/matlabcentral/answers/39522-generate-random-number-from-inverse-gamma-distribution
% rho_prop(1) = -1 + 2*rand; 
% sigma_prop(1) = 1/gamrnd(1,1/1);

rho_prop(1) = rho; 
sigma_prop(1) = sigma;
nu_prop(1) = 100;

for step = 2:steps_mcmc

        %proposed values = theta = <rho, sigma>
        rho_prop_val = -1 + 2*rand;
        sigma_prop_val = 1/gamrnd(1, 1/1);
        nu_prop_val = 1/gamrnd(1, 1/1);

        log_denominator = Ex4_Compute_Log_Likelihood(y, rho_prop(step-1), sigma_prop(step-1), beta, N, nu_prop(step-1));
        
        %1. update rho, sigma = sigma_last_prop
        log_numerator = Ex4_Compute_Log_Likelihood(y, rho_prop_val, sigma_prop(step-1), beta, N, nu_prop(step-1));
        

        if (rand <  min(1, exp(log_numerator - log_denominator)))
            rho_prop(step) = rho_prop_val;
        else
            rho_prop(step) = rho_prop(step-1);
        end

        %2. update sigma, rho = rho_last_prop
        log_numerator = Ex4_Compute_Log_Likelihood(y, rho_prop(step-1), sigma_prop_val, beta, N, nu_prop(step-1));
  
        if (rand < min(1, exp(log_numerator - log_denominator)))
            sigma_prop(step) = sigma_prop_val;
        else
            sigma_prop(step) = sigma_prop(step-1);
        end
        
        %2. update nu
        log_numerator = Ex4_Compute_Log_Likelihood(y, rho_prop(step-1), sigma_prop(step-1), beta, N, nu_prop_val);
        if (rand < min(1, exp(log_numerator - log_denominator)))
            nu_prop(step) = nu_prop_val;
        else
            nu_prop(step) = nu_prop(step-1);
        end

        if(mod(step, 10) == 0)
            fprintf('MCMC step: %i, rho = %f, sigma = %f, nu = %f \n', step, rho_prop(step-1), sigma_prop(step-1), nu_prop(step-1));
        end
end

sigma_prop = sigma_prop(sigma_prop < 10);
%MCMC_Checks(rho_prop);
%MCMC_Checks(sigma_prop);

% sigma_prop = sigma_prop(1:2500);

% rho_prop = rho_prop(1:2500);

load rho_prop_test.mat;
load sigma_prop_test.mat;
load nu_prop_test.mat;

assert(mean(sigma_prop_test == sigma_prop) == 1, 'sigma prop error');
assert(mean(rho_prop_test == rho_prop) == 1, 'rho prop error');
assert(mean(nu_prop_test == nu_prop) == 1, 'nu prop error');
