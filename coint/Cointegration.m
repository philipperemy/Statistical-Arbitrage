function [ pValue ] = Cointegration( a, b )
    ret_a = diff(log(a));
    ret_b = diff(log(b));
    
    p = polyfit(ret_b, ret_a, 1);
    beta = p(1);
    spread = a - beta * b;

    [~, pValue] = adftest(spread, 'model','TS','lags', 0);
end

% 
% length = size(a,1);
% 
% Y = [a, b];
% plot(Y);
% 
% [h,pValue,stat,cValue,reg] = egcitest(Y,'test',{'t1','t2'});
% pValue
% 
% a = reg(2).coeff(1);
% b = reg(2).coeff(2);
% plot(Y*[1;-b]-a);
% 
% [h,pValue] = adftest(Y);
% 
% p = polyfit(b,a,1);
% 
% tbl = table(a,b,'VariableNames',{'RBSA','RBSB'});
% lm = fitlm(tbl,'RBSA~RBSB-1');
% 
% ret_a = diff(log(a));
% ret_b = diff(log(b));
% tbl = table(ret_a,ret_b,'VariableNames',{'RETA','RETB'});
% lm = fitlm(tbl,'RETA~RETB-1');
% beta = lm.Coefficients.Estimate;
% spread = a - beta * b;
% 
% [h,pValue] = adftest(spread, 'model','TS','lags', 0);
% pValue
% 
% p = polyfit(ret_b, ret_a, 1);
% beta = p(1);

%500 tests per second
%780625 numbers of pairs