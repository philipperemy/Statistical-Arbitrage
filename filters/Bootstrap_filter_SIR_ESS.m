%
% Just for information. DO NOT USE IT.
% p_y_given_theta is not correct

function [estimated_states, p_y_given_theta] = Bootstrap_filter_SIR_ESS(y, rho, sigma, beta, N, steps, p_y_given_x, nu)
    
    p = zeros(N, steps);
    w = zeros(N, steps);
    p_y_given_theta = zeros(1, steps);
    estimated_states = zeros(1, steps);
    
    N_thr = round(N/3);
    
    for i = 1:steps
        %single step of sequential importance resampling
        
        if(i == 1)
            p(:,1) = randn(N,1) * sqrt(sigma^2/(1-rho^2)); % Set initial particle states 
            w(:,1) = ( 1/N );
        else
            p(:,i) = rho * p(:,i-1) + sigma * randn(N,1);
            sigma_y = beta * exp(0.5*p(:,i));
            w(:,i) = w(:,i-1) .* p_y_given_x(y(i), sigma_y, nu);
        end
        
        sum_weights = sum(w(:,i));
        p_y_given_theta(i) = (1/N)*sum_weights; %put it here or at the end?
        
        %normalize weights
        w(:,i) = w(:,i) / sum_weights;
        
        N_eff = 1/sum((w(:,i)).^2);
        if(N_eff < N_thr)
            %resampling
            %nIdx = randsample(N, N, 'true', w(:,i));
            nIdx = resampleStratified(w(:,i));
            p(:,i) = p(nIdx,i); %replace the current part set with the new one.
            w(:,i) = 1/N;
        end
        
        % Calculate state estimate
        estimated_states(i) = w(:,i)'*p(:,i);
    end
end

