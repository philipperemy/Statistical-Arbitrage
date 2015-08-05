function [  ] = ZScore( spr )
    mu = mean(spr);
    sig = std(spr);
    fprintf('mu = %f, sigma = %f\n', mu, sig);
    spr_n = (spr - mu)/sig;
    N = length(spr_n);
    plot(1:N, spr_n);
    hold on;
    plot([1 N], [-2 -2]);
    plot([1 N], [2 2]);
    hold off;
    ylim([-4 4]);
end

