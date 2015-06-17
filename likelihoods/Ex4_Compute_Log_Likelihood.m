
function [Y_log_SV, estX] = Ex4_Compute_Log_Likelihood(y, rho_prop, sigma_prop, beta, N, nu)
    steps = size(y, 2);
    [estX, ~, p_y] = Ex2_bootstrap_filter_SIR(y, rho_prop, sigma_prop, beta, N, steps, nu);
    Y_log_SV = Ex2_Estimate_Py( p_y , steps );
    Y_log_SV = Y_log_SV(length(Y_log_SV));
end