%ids = matrix of the tuples
function [ pValue, lm_formula, beta, spr, ids ] = AugmentedDickeyFullerTest( pp, id1, id2, id3 )

    s1_px = Convert(pp, id1);
    s2_px = Convert(pp, id2);
    s3_px = Convert(pp, id3);
    
    s1_name = strcat(char(pp.names(id1)), ' (', int2str(id1), ')');
    s2_name = strcat(char(pp.names(id2)), ' (', int2str(id2), ')');
    s3_name = strcat(char(pp.names(id3)), ' (', int2str(id3), ')');
    
    tbl = table(s1_px, s2_px, s3_px);
    tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3'};
    
    lm_formula = strcat(s1_name, ' ~ ', ' ', s2_name, ' + ', s3_name);
    mod = fitlm(tbl, 'stock1 ~ stock2 + stock3 - 1');
    beta = mod.Coefficients.Estimate;
    spr = s1_px - beta(1) * s2_px - beta(2) * s3_px;
    [~, pValue] = adftest(spr, 'model','TS','lags', 0);
        
    if(pValue > 0.05)

        lm_formula = strcat(s2_name, ' ~ ', s1_name, ' + ', s3_name);
        mod = fitlm(tbl, 'stock2 ~ stock1 + stock3 - 1');
        beta = mod.Coefficients.Estimate;
        spr = s2_px - beta(1) * s1_px - beta(2) * s3_px;
        [~, pValue] = adftest(spr, 'model','TS','lags', 0);

        if(pValue > 0.05)

            lm_formula = strcat(s3_name, ' ~ ', s1_name, ' + ', s2_name);
            mod = fitlm(tbl, 'stock3 ~ stock1 + stock2 - 1');
            beta = mod.Coefficients.Estimate;
            spr = s3_px - beta(1) * s1_px - beta(2) * s2_px;
            [~, pValue] = adftest(spr, 'model','TS','lags', 0);

            if(pValue > 0.05)
                lm_formula = '';
                beta = 0;
                spr = 0;
                return;
            end
        end
    end
end

