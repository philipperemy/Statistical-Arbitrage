function [ h, spread, res ] = SpreadConstructor5( stock_ids, stock_prices )

    res         = struct();
    spread      = struct();
    h           = 0;
    P_VAL_THRES = 0.05;
            
    s1_px       = stock_prices(:,1); id1 = stock_ids(1);
    s2_px       = stock_prices(:,2); id2 = stock_ids(2);
    s3_px       = stock_prices(:,3); id3 = stock_ids(3);
    s4_px       = stock_prices(:,4); id4 = stock_ids(4);
    s5_px       = stock_prices(:,5); id5 = stock_ids(5);

    %first test for a unit root. Specificity of JCI test is that order is not
    %important.
    [~,p_val_jci,~,~,~] = jcitest([s1_px, s2_px, s3_px, s4_px, s5_px], 'lags', 0);

    pval_r0 = p_val_jci.r0(1);
    pval_r1 = p_val_jci.r1(1);
    pval_r2 = p_val_jci.r2(1);
    if(pval_r0 > P_VAL_THRES)
        return;
    end
    res = struct('pval_r0_jci', pval_r0, 'pval_r1_jci', pval_r1, 'pval_r2_jci', pval_r2);

    tbl = table(s1_px, s2_px, s3_px, s4_px, s5_px);
    tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3', 'stock4', 'stock5'};
    
    %second test to determine the order of the triple. Order matters.
    mod = fitlm(tbl, 'stock1 ~ stock2 + stock3 + stock4 + stock5 - 1');
    beta = mod.Coefficients.Estimate;
    triple_output = [id1, id2, id3, id4, id5];
    spr = s1_px - (beta(1) * s2_px) - (beta(2) * s3_px) - (beta(3) * s4_px) - (beta(4) * s5_px);
    [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
        
    if(p_value_adf > P_VAL_THRES)

        mod = fitlm(tbl, 'stock2 ~ stock1 + stock3 + stock4 + stock5 - 1');
        beta = mod.Coefficients.Estimate;
        triple_output = [id2, id1, id3, id4, id5];
        spr = s2_px - (beta(1) * s1_px) - (beta(2) * s3_px) - (beta(3) * s4_px) - (beta(4) * s5_px);
        [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);

        if(p_value_adf > P_VAL_THRES)

            mod = fitlm(tbl, 'stock3 ~ stock1 + stock2 + stock4 + stock5 - 1');
            beta = mod.Coefficients.Estimate;
            triple_output = [id3, id1, id2, id4, id5];
            spr = s3_px - (beta(1) * s1_px) - (beta(2) * s2_px) - (beta(3) * s4_px) - (beta(4) * s5_px);
            [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);

            if(p_value_adf > P_VAL_THRES)
                
                mod = fitlm(tbl, 'stock4 ~ stock1 + stock2 + stock3 + stock5 - 1');
                beta = mod.Coefficients.Estimate;
                triple_output = [id4, id1, id2, id3, id5];
                spr = s4_px - (beta(1) * s1_px) - (beta(2) * s2_px) - (beta(3) * s3_px) - (beta(4) * s5_px);
                [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
                
                if(p_value_adf > P_VAL_THRES)
                   
                    mod = fitlm(tbl, 'stock5 ~ stock1 + stock2 + stock3 + stock4 - 1');
                    beta = mod.Coefficients.Estimate;
                    triple_output = [id5, id1, id2, id3, id4];
                    spr = s5_px - (beta(1) * s1_px) - (beta(2) * s2_px) - (beta(3) * s3_px) - (beta(4) * s4_px);
                    [~, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
                    
                    if(p_value_adf > P_VAL_THRES)
                        return;
                    end
                end
            end
        end
    end

    %ultimate test for a unit root
    if(~pptest(spr, 'model', 'TS', 'lags', 0))
        return;
    end
    
    %spread formed correctly
    h = 1;
    dimCol = 1;
    beta = cat(dimCol, 1, -beta); %put 1 at the beginning. -beta for the other
    spread = struct('px', spr, 'tuple', triple_output, 'beta', beta);
end
