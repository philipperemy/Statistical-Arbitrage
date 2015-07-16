%-------------------------------------------------------------------
% DO NOT USE IT PLEASE
%-------------------------------------------------------------------

function [estimated_states, estimated_states2, return_weights, p_y] = Bootstrap_filter_SIR_TwoFactors(y, rho, sigma, rho2, sigma2, beta, N, steps)
    w = zeros(N, steps);
    
    p = zeros(N, steps);
    p2 = zeros(N, steps);
    
    p_y = zeros(1, steps);
    
    return_weights = zeros(1, steps);
    
    estimated_states = zeros(1, steps);
    estimated_states2 = zeros(1, steps);
    
    FIRST_STEP = 1;
    for i=1:steps

        %draw samples from the proposal distribution (transition prior
        %probability distribution)
        if (i == FIRST_STEP)
            p(:,1) = randn(N,1) * sqrt(sigma^2/(1-rho^2)); % Set initial particle states
            p2(:,1) = randn(N,1) * sqrt(sigma2^2/(1-rho2^2)); % Set initial particle states
            w(:,1) = ( 1/N );
        else
            
            if(N_eff > 0)
                nIdx = randsample(N, N, 'true', w(:,i-1));
                p_previous = p(nIdx,i-1);
                p_previous2 = p2(nIdx,i-1);
            else
                p_previous = p(:,i-1);
                p_previous2 = p2(:,i-1);
            end
            p(:,i) = rho * p_previous + sigma * randn(N,1);
            p2(:,i) = rho2 * p_previous2 + sigma2 * randn(N,1);
        end

        %Update the weights.
        sigma_vec = beta * exp(0.5*(p(:,i)+p2(:,i)));
        w(:,i) = normpdf(y(i), 0, sigma_vec);

        return_weights(i) = sum(w(:,i));
        
        %Normalize the importance weights.
        sum_weights = sum(w(:,i));
        w(:,i) = w(:,i) / sum_weights;
        
        %Compute the effective number of particles
        N_eff = 1/sum((w(:,i)).^2);

        % Calculate state estimate
        estimated_states(i) = w(:,i)'*p(:,i);
        estimated_states2(i) = w(:,i)'*p2(:,i);
       
        % Compute p(y)
        p_y(i) = (1/N)*sum_weights;
    end

end
