function [ opt_lag ] = ADFTestOptimalLagFinder( spr )

    max_lag = round(12*(length(spr)/100)^0.25);
    stats = zeros(1, max_lag+1);
    for p = 0:max_lag
        [~,~,stat,~,~] = adftest(spr, 'model','TS','lags', p);
        stats(p+1) = stat;
    end
    
    opt_lag = find(stats == min(stats));
    fprintf('Stats min: %i\n', opt_lag);
end

