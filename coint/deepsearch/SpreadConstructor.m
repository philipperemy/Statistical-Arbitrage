
function [ res, spread, pvals_jci ] = SpreadConstructor( stock_ids, stock_prices )

    pvals_jci         = struct;
    spread            = struct;
    res               = 0; % if res = 0, spread not formed
    P_VAL_THRES       = 0.05;

    s1_px             = stock_prices(:,1); id1 = stock_ids(1);
    s2_px             = stock_prices(:,2); id2 = stock_ids(2);
    s3_px             = stock_prices(:,3); id3 = stock_ids(3);

    %first test for a unit root. Specificity of JCI test is that order is not important.
    [~,p_val_jci,~,~,~] = jcitest([s1_px, s2_px, s3_px], 'lags', 0:5);

    pval_r0 = p_val_jci.r0(1);
    pval_r1 = p_val_jci.r1(1);
    pval_r2 = p_val_jci.r2(1);
%     if(any(p_val_jci.r0 > P_VAL_THRES))
%         return;
%     end
    
    s1_px_ret = diff(log(s1_px));
    s2_px_ret = diff(log(s2_px));
    s3_px_ret = diff(log(s3_px));
    
    tbl = table(s1_px_ret, s2_px_ret, s3_px_ret);
    tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3'};
    
    %second test to determine the order of the triple. Order matters.
    mod = fitlm(tbl, 'stock1 ~ stock2 + stock3 - 1');
    beta = mod.Coefficients.Estimate;
    triple_output = [id1, id2, id3];
    spr = s1_px - beta(1) * s2_px - beta(2) * s3_px;
    [~, p_value_adf] = adftest(spr, 'model','ARD','lags', 0);
        
    if(p_value_adf > P_VAL_THRES)

        mod = fitlm(tbl, 'stock2 ~ stock1 + stock3 - 1');
        beta = mod.Coefficients.Estimate;
        triple_output = [id2, id1, id3];
        spr = s2_px - beta(1) * s1_px - beta(2) * s3_px;
        [~, p_value_adf] = adftest(spr, 'model','ARD','lags', 0);

        if(p_value_adf > P_VAL_THRES)

            mod = fitlm(tbl, 'stock3 ~ stock1 + stock2 - 1');
            beta = mod.Coefficients.Estimate;
            triple_output = [id3, id1, id2];
            spr = s3_px - beta(1) * s1_px - beta(2) * s2_px;
            [~, p_value_adf] = adftest(spr, 'model','ARD','lags', 0);

            if(p_value_adf > P_VAL_THRES)
                return;
            end
        end
    end

    %ultimate test for a unit root
    if(~pptest(spr, 'model', 'ARD', 'lags', 0))
        return;
    end
    
%     if(~vratiotest(spr))
%        return; 
%     end
    
    %spread formed correctly
    res = 1;
    dimCol = 1;
    beta = cat(dimCol, 1, -beta); %put 1 at the beginning. -beta for the other
    spread = struct('px', spr, 'tuple', triple_output, 'beta', beta); 
    pvals_jci = struct('pval_r0_jci', pval_r0, 'pval_r1_jci', pval_r1, 'pval_r2_jci', pval_r2);
end
