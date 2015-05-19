
function [Y_log_SV] = Compute_Log_Likelihood_NormalLeverage(y, rho_prop, sigma_prop, beta, N, cor)
    steps = size(y, 2);
    [~, ~, p_y, error] = Bootstrap_filter_SIR_NormalLeverage(y, rho_prop, sigma_prop, beta, N, steps, cor);
    
    if(error == 1)
        Y_log_SV = -Inf;
        return;
    end
    
    Y_log_SV = Ex2_Estimate_Py( p_y , steps );
    Y_log_SV = Y_log_SV(length(Y_log_SV));
end