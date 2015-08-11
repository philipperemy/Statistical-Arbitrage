function [ ] = Draw_Bollinger_Bands( range, wsize, nstd, dis_HIGH, dis_LOW, dis_CLOSE, dis_OPEN, hasLegend )
    figure;
    [mid, uppr, lowr] = bollinger(dis_CLOSE(range), wsize, 0, nstd);
    plot([Fit_NA(mid) Fit_NA(uppr) Fit_NA(lowr)], 'LineWidth', 2);
    hold on;
    candle(dis_HIGH(range), dis_LOW(range), dis_CLOSE(range), dis_OPEN(range), 'black');
    if(hasLegend)
        legend('mid band', 'upper band', 'lower band', 'stock price');
    end
    hold on;
    xlim([wsize range(end)-range(1)]);
    ylabel('Prices (USD)');
    xlabel('Time (Days)');

end

