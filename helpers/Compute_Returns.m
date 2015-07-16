function [ returns ] = Compute_Returns( prices )
    N = length(prices);
    returns = zeros(1, N);
    returns(1) = 0;
    for i = 2:N
        returns(i) = prices(i)/prices(i-1)-1;
    end
end

