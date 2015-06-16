function [ p_y ] = GelfandDey_LogMarginalY(st, thetas )

    N = size(thetas,1);
    log_ratio = zeros(1, N);
    for i = 1:N
        rho = thetas(i, 1);
        sigma = thetas(i, 2);
        nu = thetas(i, 3);
        
        log_g_theta = log(mvnpdf([rho sigma nu], [rho sigma nu], eye(3))); %order of parameter is not important
        log_p_theta = log((1/2)*gampdf(sigma,1,1)*gampdf(nu,1,1));
        log_p_y_theta = Ex4_Compute_Log_Likelihood(st.y, rho, sigma, 0.5, 1000, nu);
        
        log_ratio(i) = log_g_theta - log_p_theta - log_p_y_theta;
    end
    
    max_log_ratio = max(log_ratio);
    log_ratio_centered = log_ratio - max_log_ratio;
    p_y = max_log_ratio + log(sum(exp(log_ratio_centered)));
    p_y = - p_y;
end

