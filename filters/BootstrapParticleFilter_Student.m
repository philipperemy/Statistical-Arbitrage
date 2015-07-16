
function [log_p_y_given_theta, estimated_states] = BootstrapParticleFilter_Student(y, rho, sigma, beta, nu, N, p_y_given_x)
    
    %fprintf('[PF] rho = %f, sigma = %f, beta = %f, nu = %f\n', rho, sigma, beta, nu);
    T = length(y);
    p = zeros(N, T);
    w = zeros(N, T);
    estimated_states = zeros(1,T);
    log_p_y_given_theta = zeros(1,T);
    
    p(:,1)              =   randn(N,1) * sqrt(sigma^2/(1-rho^2));
    w(:,1)              =   p_y_given_x(y(1), beta*exp(0.5*p(:,1)), nu);
    log_p_y_given_theta(1)  =   log(mean(w(:,1)));
    
    for t = 2:T

        try
            nIdx = randsample(N, N, 'true', w(:,t-1));
        catch
            nIdx = 1:N;
        end
        
        p(:,t)             = rho * p(nIdx,t-1) + sigma * randn(N,1);
        w(:,t)             = p_y_given_x(y(t), beta * exp(0.5*p(:,t)), nu);
        log_p_y_given_theta(t) = log_p_y_given_theta(t-1) + log(mean(w(:,t)));
        
        % Calculate state estimate
        norm_w = w(:,t) / sum(w(:,t));
        estimated_states(t) = norm_w'*p(:,t);
    end
    
    log_p_y_given_theta = log_p_y_given_theta(end);
end