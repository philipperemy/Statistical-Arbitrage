function [ normalized_std_diff ] = Volatility_Estimation_Bands_Mean_Zero( spread, returns_volatility_var)
    prices = spread.px;
    returns_volatility_sd = sqrt(returns_volatility_var);
    M = 1000; %generate M observations of the process
    T = length(prices);
    generated_processes = zeros(M, T);
    vol_returns = zeros(T,1);
    for t = 2:T
        vol_returns(t) = returns_volatility_sd(t)*prices(t-1);
        generated_processes(:,t) = normrnd(prices(t-1), vol_returns(t), M, 1);
        %generated_processes(:,t) = prices(t-1)*(1+returns_volatility_sd(t)*randn(N,1));
    end

    beg = 4;
    figure; hold on;
    plot(beg:length(prices), prices(beg:end)', 'Color', 'r', 'LineWidth', 1);
    plot(beg:length(prices), generated_processes(:, beg:end)', 'Color',[0,0.7,0.9]);
    legend('Generated prices trajectories from SV', 'Observed spread prices');
    ylim([min(prices)-5 max(prices)+5]);
    xlabel('Time (days)');
    ylabel('Price (USD)');
    plot(beg:length(prices), prices(beg:end)', 'Color', 'r', 'LineWidth', 1);
    hold off;
    
    figure;
    plot(vol_returns(beg:end).^2);
    xlabel('Time (days)');
    ylabel('Price (USD)');
    legend('Conditional variance Var(S_t | S_{t-1}, \theta)');

    mvstd_period = 20;
    movstd_mat = zeros(M,T);
    for n = 1:M
        movstd_mat(n,:) = movingstd(generated_processes(n,:), mvstd_period);
    end

    eval_vol    = mean(movstd_mat, 1)';
    common_vol  = movingstd(prices, mvstd_period);
    Y = [eval_vol common_vol];
    
    %http://doc.instantreality.org/tools/color_calculator/
    figure; hold on;
    Y = Y(mvstd_period+5:end,:);
    plot(Y(:,1), 'LineWidth', 1, 'Color', [1 .5 0]); %trick for legend
    plot(Y(:,2), 'LineWidth', 1, 'Color', [0.956, 0.886, 0.243]);
    plot(movstd_mat(:,mvstd_period+5:end)',  'Color',[0.4,0.4,0.4]);
    legend('Aggregated moving std from SV', 'Moving std from standard', 'Monte Carlo moving std from SV');
    plot(Y(:,1), 'LineWidth', 1, 'Color', [1 .5 0]);
    plot(Y(:,2), 'LineWidth', 1, 'Color', [0.956, 0.886, 0.243]);
    xlabel('Time (days)');
    ylabel('Price (USD)');
    hold off;
    
    figure; hold on;
    plot(Y(:,2), 'LineWidth', 1, 'Color', [0,0.7,0.9]);
    plot(quantile(movstd_mat(:,mvstd_period+5:end), 0.95)', '--', 'LineWidth', 2);
    plot(quantile(movstd_mat(:,mvstd_period+5:end), 0.05)', '--', 'LineWidth', 2);
    legend('Moving std from standard', 'Monte Carlo Quantile 0.95', 'Monte Carlo Quantile 0.05');
    xlabel('Time (days)');
    ylabel('Price (USD)');
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bollinger Bands
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ma = tsmovavg(prices,'s',mvstd_period,1);
    d = 1.5;
    Y = [ma-d*eval_vol ma prices ma+d*eval_vol];
    Y = [Y ma-d*common_vol ma+d*common_vol]; %concat

    figure;
    plot(Y(mvstd_period+5:end,:));
    legend('boll-e', 'ma', 'prices', 'boll+e', 'boll-r', 'boll+r');
    xlabel('Time (days)');
    ylabel('Price (USD)');
    normalized_std_diff = (eval_vol-common_vol)./prices;
    figure;
    plot(normalized_std_diff(30:end));

end

