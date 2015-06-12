function [ profits, vols, sharpes, I ] = CrossValidationParallel( pp, Spread, strat_func, nstd_range, wts_range, wsize_range )
    %tic;
    I = combvec(nstd_range, wts_range, wsize_range)';
    profits = zeros(1, length(I));
    vols = zeros(1, length(I));
    sharpes = zeros(1, length(I));
    parfor j = 1:length(I);
        tmp = I(j,:);
        wsize = tmp(3);
        wts = tmp(2);
        nstd = tmp(1);
        boll_conf = struct('wsize', wsize, 'wts', wts, 'nstd', nstd);
        [pl, balance_cum] = strat_func( pp, Spread, 1, length(Spread.px), boll_conf, 0 );
        sharpes(j) = (mean(pl)/std(pl))*sqrt(252); % 0.7776 http://www.nuclearphynance.com/Show%20Post.aspx?PostIDKey=167056
        profits(j) = balance_cum(end);
        vols(j) = std(pl);
    end
    %toc;
    
end

