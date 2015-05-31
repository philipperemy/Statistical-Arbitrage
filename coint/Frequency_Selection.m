%ids = matrix of the tuples
function [  ] = Frequency_Selection( pp, ids, ma_period )
    for i = 1:size(ids,1)
        s1 = ids(i, 1);
        s2 = ids(i, 2);
        s3 = ids(i, 3);
        
        s1_px = Convert(pp, s1);
        s2_px = Convert(pp, s2);
        s3_px = Convert(pp, s3);
        
        tbl = table(s1_px, s2_px, s3_px);
        tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3'};
        
        %
        % [h,pValue,stat,cValue,reg] = pptest(A, 'model', 'TS')
        % pvalue less than 0.05. Trend stationary
        % http://fr.mathworks.com/help/econ/pptest.htmls
        % RW : very high pvalue
        %
        
        %Optimal lag
        %http://fr.mathworks.com/help/econ/var-models.html#bswxr8u-27
        
        
        case_id = 1;
        lm_formula = strcat(char(pp.names(s1)), ' ~ ', ' ', char(pp.names(s2)), ' + ', char(pp.names(s3)));
        mod = fitlm(tbl, 'stock1 ~ stock2 + stock3 - 1');
        beta = mod.Coefficients.Estimate;
        spr = s1_px - beta(1) * s2_px - beta(2) * s3_px;
        [~, pValue] = adftest(spr, 'model','TS','lags', 0);
        
        if(pValue > 0.05)
            
            lm_formula = strcat(char(pp.names(s2)), ' ~ ', char(pp.names(s1)), ' + ', char(pp.names(s3)));
            mod = fitlm(tbl, 'stock2 ~ stock1 + stock3 - 1');
            beta = mod.Coefficients.Estimate;
            spr = s2_px - beta(1) * s1_px - beta(2) * s3_px;
            
            [~, pValue] = adftest(spr, 'model','TS','lags', 0);
            case_id = 2;
            if(pValue > 0.05)
                case_id = 3;
                
                lm_formula = strcat(char(pp.names(s3)), ' ~ ', char(pp.names(s1)), ' + ', char(pp.names(s2)));
                mod = fitlm(tbl, 'stock3 ~ stock1 + stock2 - 1');
                beta = mod.Coefficients.Estimate;
                spr = s3_px - beta(1) * s1_px - beta(2) * s2_px;
                
                [~, pValue] = adftest(spr, 'model','TS','lags', 0);
                if(pValue > 0.05)
                    warning('error should not happen\n');
                    continue;
                end
            end
        end
        
        freq = Revertion_Frequency(spr, ma_period);
        fprintf('%i, %s, freq = %f, reversion period = %i \n', case_id,...
            lm_formula, freq, round(1/freq));
    end
end

