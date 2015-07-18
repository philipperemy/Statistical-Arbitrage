function [ log_marginal_likelihood_mat ] = ParticleCountProfilingStudentSV_Unit( st, N_seq, MC_count, rho, sigma, beta, nu)
    j = 1;
    p_y_given_x = @(y,sig, nu) tpdf(y./sig, nu)./sig;
    log_marginal_likelihood_mat = zeros(length(N_seq), MC_count);
    for N = N_seq
        for i = 1:MC_count
            log_marginal_likelihood_mat(j, i) = BootstrapParticleFilter_Student(st.y, rho, sigma, beta, nu, N, p_y_given_x);
            fprintf('%i, %i\n', i, N);
        end
        j = j + 1;
    end
end