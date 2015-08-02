function [ res ] = ADFTestOptimalLagFinder(pp, id1, id2, id3, date_range)

    Y = [GetPrice(pp, id1), GetPrice(pp, id2), GetPrice(pp, id3)];
    if(exist('range', 'var'))
        Y = Y(date_range,:);
    end

    [res, spread, ~] = SpreadConstructor( [id1 id2 id3], Y);
    if (res == 0)
        res = -1;
        return;
    end
    spr = spread.px;
    
    max_lag = round(12*(length(Y)/100)^0.25);
    for p = 0:max_lag
        [~,~,stat,~,reg] = adftest(spr, 'model','TS','lags', p);
        aic(p+1)   = reg.AIC;
        bic(p+1)   = reg.BIC;
        ll(p+1)    = reg.LL;
        hqc(p+1)   = reg.HQC;
        stats(p+1) = stat;
        fprintf('computing for lag %i\n', p);
    end
    
    res = find(stats == min(stats));
    fprintf('Stats min: %i\n', res);
    fprintf('BIC min: %i\n', find(bic == min(bic)));
    fprintf('LL min: %i\n', find(ll == min(ll)));
    fprintf('HQC min: %i\n', find(hqc == min(hqc)));
    fprintf('AIC min: %i\n', find(aic == min(aic)));

    subplot(2,1,1);
    plot(stats);
    title('t-stat ADF');
    subplot(2,1,2);
    plot(spr);
    title('spread');
    
%     subplot(3,2,1);
%     plot(aic);
%     title('aic');
%     subplot(3,2,2);
%     plot(bic);
%     title('bic');
%     subplot(3,2,3); 
%     plot(ll);
%     title('ll');
%     subplot(3,2,4);
%     plot(hqc);
%     title('hqc');
%     subplot(3,2,5);
%     plot(stats);
%     title('stats');
%     subplot(3,2,6);
%     plot(spr);
%     title('spread');

end

