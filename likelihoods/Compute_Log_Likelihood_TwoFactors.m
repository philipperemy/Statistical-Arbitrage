
function [Y_log_SV] = Compute_Log_Likelihood_TwoFactors(y, rho_prop, sigma_prop, rho_prop2, sigma_prop2, beta, N)
    steps = size(y, 2);
    [~, ~, ~, p_y] = Bootstrap_filter_SIR_TwoFactors(y, rho_prop, sigma_prop, rho_prop2, sigma_prop2, beta, N, steps);
    Y_log_SV = Ex2_Estimate_Py( p_y , steps );
    Y_log_SV = Y_log_SV(length(Y_log_SV));
end