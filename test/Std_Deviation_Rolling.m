function [ out ] = Std_Deviation_Rolling( prices, period )
    T = length(prices);
    out = zeros(1, T);
    for i = (period+1):T
       out(i) = std(prices((i-period):i)); 
    end
end

% for idx = wsize:length(nnandata)
%    mstd(idx-wsize+1, :) = std(nnandata(idx-wsize+1:idx));
% end