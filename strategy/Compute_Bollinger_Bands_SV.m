function [ mid, uppr, lowr ] = Compute_Bollinger_Bands_SV(data, wsize, wts, nstd, returns_volatility_var)
    %wts is not used. Not yet.
    returns_volatility_sd = sqrt(returns_volatility_var);
    M = 100;%for speed. MC param
    T = length(data);
    generated_processes = zeros(M, T);
    vol_returns = zeros(T,1);
    for t = 2:T
        vol_returns(t) = returns_volatility_sd(t)*data(t-1);
        generated_processes(:,t) = normrnd(data(t-1), vol_returns(t), M, 1);
    end
    
    movstd_mat = zeros(M,T);
    for n = 1:M
        movstd_mat(n,:) = movingstd(generated_processes(n,:), wsize, 'backward');
    end

    sv_vol = mean(movstd_mat, 1)';
    
    if(wts == 1)
        [mid,~,~] = bollinger(data, wsize, wts, nstd);
    else
         mid = tsmovavg(data,'s',wsize,1);
    end
    
    uppr   = mid+nstd*sv_vol;
    lowr   = mid-nstd*sv_vol;
end