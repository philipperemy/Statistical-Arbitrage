function [ ] = MCMC_Checks( props )
    
    subplot(2,2,1);
    plot(props);
    ylim([0 max(props)]);
    
    subplot(2,2,2);
    qqplot(props);
    
    subplot(2,2,3); 
    hist(props, 50);
    
    subplot(2,2,4);
    autocorr(props, 100);
    
    acc = mean(diff(props,1) ~= 0);
    fprintf('Mean : %f, Median : %f, Min : %f, Max : %f, Acceptance rate ratio : %f \n', mean(props), median(props), min(props), max(props), acc);
end

