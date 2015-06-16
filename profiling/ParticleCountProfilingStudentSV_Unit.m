function [ N_ideal ] = ParticleCountProfilingStudentSV_Unit( st, N_seq, MC_count, rho, sigma, beta, nu, it)
    %N_seq = 10:10:N_max;
    j = 1;
    log_marginal_likelihood_mat = zeros(length(N_seq), MC_count);
    for N = N_seq
        for i = 1:MC_count
            log_marginal_likelihood_mat(j, i) = Ex4_Compute_Log_Likelihood(st.y, rho, sigma, beta, N, nu);
            fprintf('%i, %i, %i\n', i, N, it);
        end
        j = j + 1;
    end
    true_val = mean(log_marginal_likelihood_mat, 2);
    true_val = true_val(end);
    
    threshold = 0.5*abs(true_val)/st.steps;
    ret = std(log_marginal_likelihood_mat,0,2);
    N_ideal = find(ret < threshold == 1);
    N_ideal = N_ideal(1);
    
    if(numel(N_ideal) == 0)
        N_ideal = N_seq(end);
    end
    
    N_ideal = N_seq(N_ideal);
end