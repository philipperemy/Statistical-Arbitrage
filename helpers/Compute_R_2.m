function [ val ] = Compute_R_2( stock_1, stock_2, stock_3, stock_4 )
    A = [stock_1 stock_2 stock_3 stock_4];
    cor_mat = corr(A);
    c = cor_mat(2:4,1);
    A(:,1) = [];
    R = corr(A);
    val = c'*inv(R)*c;
end

%   LM = FITLM(X,Y) fits a linear regression model using the column vector
%   Y as a response variable and the columns of the matrix X as predictor
%   variables, and returns the linear model LM.
%   mdl = fitlm([stock_2 stock_3 stock_4], stock_1, 'Intercept', false);
%   val = mdl.Rsquared.Ordinary;