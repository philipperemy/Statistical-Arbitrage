
%Works now for lags = 0 ONLY
function [ Spreads, beta, tuples_ret ] = SpreadConstructor( pp, tuples )
    c = 1;
    for i = 1:size(tuples,1)
        s1 = tuples(i, 1);
        s2 = tuples(i, 2);
        s3 = tuples(i, 3);
        
        s1_px = Convert(pp, s1);
        s2_px = Convert(pp, s2);
        s3_px = Convert(pp, s3);
        
        tbl = table(s1_px, s2_px, s3_px);
        tbl.Properties.VariableNames = { 'stock1', 'stock2', 'stock3'};
        
        [ pValue, lm_formula, beta, spr, tuples_ret ]  = AugmentedDickeyFullerTest(pp, s1, s2, s3);

        if(pValue > 0.05)
            continue;
        end
        
        period = 20;
        freq = Revertion_Frequency( spr, period );
        Spreads(c) = struct('px', spr, 'formula', lm_formula, 'beta', beta, 'freq', freq,...
            'period', period);
        fprintf('Spread [%s] Constructed\n', lm_formula);
        c = c + 1;
    end
end

%
% [h,pValue,stat,cValue,reg] = pptest(A, 'model', 'TS')
% pvalue less than 0.05. Trend stationary
% http://fr.mathworks.com/help/econ/pptest.htmls
% RW : very high pvalue
%

%Optimal lag
%http://fr.mathworks.com/help/econ/var-models.html#bswxr8u-27
%http://stats.stackexchange.com/questions/21539/what-is-the-correct-procedure-to-choose-the-lag-when-performing-johansen-cointeg