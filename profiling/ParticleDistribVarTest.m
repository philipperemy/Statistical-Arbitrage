function [ var_part_dist, N_seq ] = ParticleDistribVarTest( steps )
    MC_steps = 10;
    N_seq = 100:100:3000;
    rho = 0.91;       %Scale parameter for X
    sigma = 1;        %Standard deviation of the X (process)
    beta = 0.5;       %Scale parameter for Y
    nu = 3;           %Degrees of freedom of the innovation of the measurement
    st = StochasticVolatilityModelStudent(rho, sigma, beta, steps, nu);
    particles_dist_test = ParticleCountProfilingStudentSV_Unit( st, N_seq, MC_steps, rho, sigma, beta, nu );
    var_part_dist = var(particles_dist_test, 0, 2);
end