function [ log_marginal_likelihood_mat ] = Sensibility_with_rho( N_opt, rho_seq, rho, sigma, beta, nu, steps )
    st = StochasticVolatilityModelStudent(rho, sigma, beta, steps, nu);
    MAX = 100;
    j = 1;
    log_marginal_likelihood_mat = zeros(length(rho_seq), MAX);
    for rho = rho_seq
        log_marginal_likelihood_mat(j, :) = ParticleCountProfilingStudentSV_Unit( st, N_opt, MAX, rho, sigma, beta, nu );
        fprintf('%rho = %f\n', rho);
        j = j + 1;
    end
end

