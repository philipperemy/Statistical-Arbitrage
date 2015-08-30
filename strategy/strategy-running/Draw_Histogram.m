function [ ] = Draw_Histogram( vals, xlimit, bins, legend_name, legend_location )
    figure;
    hist(vals, bins);
    mu = mean(vals);
    if(exist('xlimit', 'var'))
        xlim(xlimit);
    end

    hold on;
    plot([mu,mu],ylim,'r--','LineWidth',2);
    
    if(exist('legend_name', 'var'))
       legend(legend_name, 'mean','Location', legend_location);
    end
    xlabel('Value');
    ylabel('Frequency');
    hold off;
    drawnow;
end

