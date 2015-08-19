
function [log_p_y_given_theta, estimated_states, estimated_states2] = BootstrapParticleFilter_TwoFactorsCor(y, rho1, sigma1, rho2, sigma2, beta, cor, N, p_y_given_x)
    
    %fprintf('rho1 = %f, sigma1 = %f, rho2 = %f, sigma2 = %f, beta = %f, cor = %f\n', rho1, sigma1, rho2, sigma2, beta, cor);
    T                       = length(y);
    p1                      = zeros(N, T);
    p2                      = zeros(N, T);
    w                       = zeros(N, T);
    estimated_states        = zeros(1,T);
    estimated_states2       = zeros(1,T);
    log_p_y_given_theta     = zeros(1,T);
    
    p1(:,1)                 =   randn(N,1) * sqrt(sigma1^2/(1-rho1^2));
    p2(:,1)                 =   randn(N,1) * sqrt(sigma2^2/(1-rho2^2)); % Set initial particle states
    w(:,1)                  =   p_y_given_x(y(1), beta * exp( 0.5*( p1(:,1) + p2(:,1) )) );
    log_p_y_given_theta(1)  =   log(mean(w(:,1)));
    norm_w                  =   (w(:,1) / sum(w(:,1)));
    estimated_states(1)     =   norm_w'*p1(:,1);
    estimated_states2(1)    =   norm_w'*p2(:,1);
	last_innov_y            =   y(1) / ( beta * exp(0.5*estimated_states(1) + 0.5*estimated_states2(1)) );
    
    for t = 2:T

        try
            nIdx = randsample(N, N, 'true', w(:,t-1));
        catch
            nIdx = 1:N;
        end
        
		innov_x                 = sigma1 * (cor * last_innov_y + sqrt(1-cor^2)*randn(N,1));
        p1(:,t)                 = rho1 * p1(nIdx,t-1) + innov_x;
        p2(:,t)                 = rho2 * p2(nIdx,t-1) + sigma2 * randn(N,1);
        
        w(:,t)                  = p_y_given_x(y(t), beta * exp( 0.5*( p1(:,t) + p2(:,t) )));
        log_p_y_given_theta(t)  = log_p_y_given_theta(t-1) + log(mean(w(:,t)));
        
        % Calculate state estimate
        norm_w                  = w(:,t) / sum(w(:,t));
        estimated_states(t)     = norm_w'*p1(:,t);
        estimated_states2(t)    = norm_w'*p2(:,t);
		last_innov_y            = y(t) / ( beta * exp(0.5*estimated_states(t) + 0.5*estimated_states2(t)) );
    end
    
    log_p_y_given_theta = log_p_y_given_theta(end);
end