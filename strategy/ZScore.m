function [  ] = ZScore( spr, zscore_conf )
    mu = mean(spr);
    sig = std(spr);
    fprintf('mu = %f, sigma = %f\n', mu, sig);
    spr_n = (spr - mu)/sig;
    N = length(spr_n);
    plot_line_width = 1;
    subplot(2,1,1);
    plot(1:N, spr, 'Color', 'Blue');
    ylabel('S(t)');
    xlabel('Time (Days)');
    subplot(2,1,2);
    hold on;
    plot(1:N, spr_n, 'Color', 'Blue');
    plot([1 N], [zscore_conf.sell_open zscore_conf.sell_open], 'LineWidth', plot_line_width);
    plot([1 N], [zscore_conf.sell_close zscore_conf.sell_close], 'LineWidth', plot_line_width);
    plot([1 N], [zscore_conf.buy_open zscore_conf.buy_open], 'LineWidth', plot_line_width);
    plot([1 N], [zscore_conf.buy_close zscore_conf.buy_close], 'LineWidth', plot_line_width);
    ylabel('Z(t)');
    xlabel('Time (Days)');
    ylim([-4 4]);
    plot([1 N], [4 4], 'Color', 'black');
    hold off;
end

