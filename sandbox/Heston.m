clear;
spot = 50;
v0 = 1;

max1 = 100; 
alpha = 0.01; 
beta = 0.001; 
sigma = 0.075; 
theta = 0.001; 
rho = 0.04; 
t = 90;
r = 0.02;
s(1) = spot; 
v(1) = v0; 
M = [1 rho; rho 1]; % The diagonal matrix for use in cholesky method 

for i=2:(t+2) 
    epsilon = chol(M)'*randn(2,1); 
    % Generates random numbers in order to determine epsilons by use of 
    %cholesky method 
    deltatime =1/252; 
    % time discretization over the life of the option' the time step length is 
    % one day. 
    v(i) = v(i-1) + ( alpha - beta * v(i-1)) * deltatime + sigma * ((v(i-1) * deltatime)^.5) * epsilon(2); 
    s(i) = s(i-1) + (r- theta) * s(i-1)      * deltatime + s(i-1) * ((v(i-1) * deltatime)^.5) * epsilon(1); 
end 

plot(v)
figure;
plot(s)
