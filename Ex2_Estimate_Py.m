function [ Y_log_SV ] = Ex2_Estimate_Py( estimate_py, steps )
    Y_log_SV = zeros(1, steps);
    for i = 2:steps
        Y_log_SV(i) = sum(log(estimate_py(2:i)));
    end
end

