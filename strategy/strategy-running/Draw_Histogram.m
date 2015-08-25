function [ ] = Draw_Histogram( vals, xlimit )
    figure;
    hist(vals, 40);
    mu = mean(vals);
    xlim(xlimit);
    hold on;
    plot([mu,mu],ylim,'r--','LineWidth',2);
    hold off;
    drawnow;
end

