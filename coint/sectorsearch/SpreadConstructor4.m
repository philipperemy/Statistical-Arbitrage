function [ h, spread, res ] = SpreadConstructor4( stock_ids, stock_prices, beta )

    res         = struct();
    spread      = struct();
    h           = 0;
    P_VAL_THRES = 0.05;
            
    s1_px       = stock_prices(:,1);
    s2_px       = stock_prices(:,2);
    s3_px       = stock_prices(:,3);
    s4_px       = stock_prices(:,4);
    
    %first test each one for unit root
    h1 = adftest(s1_px, 'model', 'TS', 'lags', 0:2);
    h2 = adftest(s2_px, 'model', 'TS', 'lags', 0:2);
    h3 = adftest(s3_px, 'model', 'TS', 'lags', 0:2);
    h4 = adftest(s4_px, 'model', 'TS', 'lags', 0:2);
    
    if(any(h1) || any(h2) || any(h3) || any(h4)) %at least one is 1. which means TS.
       return; 
    end
    
    if(exist('beta', 'var'))
        h = 1;
        spr = (beta(1) * s1_px) + (beta(2) * s2_px) + (beta(3) * s3_px) + (beta(4) * s4_px);
        spread = struct('px', spr, 'tuple', stock_ids, 'beta', beta);
        return;
    end

    %first test for a unit root. Specificity of JCI test is that order is not
    %important.
    [~,p_val_jci,~,~,~] = jcitest([s1_px, s2_px, s3_px, s4_px], 'lags', 0);

    pval_r0 = p_val_jci.r0(1);
    pval_r1 = p_val_jci.r1(1);
    pval_r2 = p_val_jci.r2(1);
    if(pval_r0 > P_VAL_THRES)
        return;
    end
    res = struct('pval_r0_jci', pval_r0, 'pval_r1_jci', pval_r1, 'pval_r2_jci', pval_r2);

    tbl = table(s1_px, s2_px, s3_px, s4_px);
    tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3', 'stock4'};
    
    %second test to determine the order of the triple. Order matters.
    mod = fitlm(tbl, 'stock1 ~ stock2 + stock3 + stock4 - 1');
    beta = mod.Coefficients.Estimate;
    spr = s1_px - (beta(1) * s2_px) - (beta(2) * s3_px) - (beta(3) * s4_px);
    [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
        
    if(p_value_adf > P_VAL_THRES)
        return;
    end

    %ultimate test for a unit root
    if(~pptest(spr, 'model', 'TS', 'lags', 0))
        return;
    end
    
    %spread formed correctly
    h = 1;
    dimCol = 1;
    beta = cat(dimCol, 1, -beta); %put 1 at the beginning. -beta for the other
    spread = struct('px', spr, 'tuple', stock_ids, 'beta', beta);
end
