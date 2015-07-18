
function [log_p_y_given_theta, estimated_states] = BootstrapParticleFilter_NormalLeverage(y, rho, sigma, beta, cor, mu, N, p_y_given_x)
    
    %fprintf('[PF] rho = %f, sigma = %f, beta = %f, cor = %f, mu = %f\n', rho, sigma, beta, cor, mu);
    T = length(y);
    p = zeros(N, T);
    w = zeros(N, T);
    estimated_states = zeros(1,T);
    log_p_y_given_theta = zeros(1,T);
    
    p(:,1)                  =   randn(N,1) * sqrt(sigma^2/(1-rho^2));
    w(:,1)                  =   p_y_given_x(y(1), beta*exp(0.5*p(:,1)));
    log_p_y_given_theta(1)  =   log(mean(w(:,1)));
    estimated_states(1)     =   (w(:,1) / sum(w(:,1)))'*p(:,1); %correct it in other filters.
    last_innov_y            =   y(1) / ( beta * exp(estimated_states(1)/2) );
    
    for t = 2:T

        try
            nIdx = randsample(N, N, 'true', w(:,t-1));
        catch
            nIdx = 1:N;
        end
        %this.mu + this.rho  * (x(i)-this.mu)
        innov_x                 = sigma * (cor * last_innov_y + sqrt(1-cor^2)*randn(N,1));
        p(:,t)                  = mu + rho * (p(nIdx,t-1)-mu) + innov_x;
        w(:,t)                  = p_y_given_x(y(t), beta * exp(0.5*p(:,t)));
        log_p_y_given_theta(t)  = log_p_y_given_theta(t-1) + log(mean(w(:,t)));
        
        % Calculate state estimate
        norm_w                  = w(:,t) / sum(w(:,t));
        estimated_states(t)     = norm_w'*p(:,t);
        last_innov_y            = y(t) / ( beta * exp(estimated_states(t)/2) );
       
    end
    
    log_p_y_given_theta = log_p_y_given_theta(end);
end