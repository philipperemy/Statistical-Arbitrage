function [  ] = ZScore( spr, zscore_conf )
    mu = mean(spr);
    sig = std(spr);
    fprintf('mu = %f, sigma = %f\n', mu, sig);
    spr_n = (spr - mu)/sig;
    N = length(spr_n);
    plot(1:N, spr_n);
    hold on;
    plot([1 N], [zscore_conf.sell_open zscore_conf.sell_open]);
    plot([1 N], [zscore_conf.sell_close zscore_conf.sell_close]);
    plot([1 N], [zscore_conf.buy_open zscore_conf.buy_open]);
    plot([1 N], [zscore_conf.buy_close zscore_conf.buy_close]);
    hold off;
    ylim([-4 4]);
 
end

