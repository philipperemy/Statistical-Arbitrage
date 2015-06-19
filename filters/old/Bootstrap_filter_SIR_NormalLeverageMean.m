%-------------------------------------------------------------------
% Bootstrap filter Exercise 2
%-------------------------------------------------------------------

function [estimated_states, return_weights, p_y, error] = Bootstrap_filter_SIR_NormalLeverageMean(y, rho, sigma, beta, N, steps, cor, mu)
    w = zeros(N, steps);
    p = zeros(N, steps);
    p_y = zeros(1, steps);
    return_weights = zeros(1, steps);
    estimated_states = zeros(1, steps);
    FIRST_STEP = 1;
    innov_Y = 1;
    error = 0;
    for i=1:steps

        %draw samples from the proposal distribution (transition prior
        %probability distribution)
        if (i == FIRST_STEP)
            p(:,1) = randn(N,1) * sqrt(sigma^2/(1-rho^2)); % Set initial particle states
            w(:,1) = ( 1/N );
        else
            
            if(N_eff > 0)
            %if(N_eff < 500) %arbitrary value => if you use it. Make sure
            %you reweights to 1/N. and use them in the update weights by
            %exp( log(w:,i-1) + ...)
                %Resampling phase
                nIdx = randsample(N, N, 'true', w(:,i-1));
                p_previous = p(nIdx,i-1);
            else
                p_previous = p(:,i-1);
            end
            p(:,i) = mu + rho * (p_previous-mu) + sigma * (cor * innov_Y + sqrt(1-cor^2)*randn(N,1));
            
        end

        %Update the weights.
        sigma_vec = beta * exp(0.5*p(:,i));
        w(:,i) = normpdf(y(i), 0, sigma_vec);
        
        return_weights(i) = sum(w(:,i));
        
        %Normalize the importance weights.
        sum_weights = sum(w(:,i));
        
        if(sum_weights == 0)
            error = 1;
            return;
        end
        
        w(:,i) = w(:,i) / sum_weights;
        
        %Compute the effective number of particles
        N_eff = 1/sum((w(:,i)).^2);

        % Calculate state estimate
        estimated_states(i) = w(:,i)'*p(:,i);
        
        innov_Y = y(i) / ( beta * exp(estimated_states(i)/2) );
       
        % Compute p(y)
        p_y(i) = (1/N)*sum_weights;
    end

end
