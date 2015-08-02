function [ val ] = Compute_R_2_ETF( stock_etf, stock_1, stock_2, stock_3, stock_4 )
    A = [stock_etf stock_1 stock_2 stock_3 stock_4];
    cor_mat = corr(A);
    c = cor_mat(2:5,1);
    A(:,1) = []; %remove first column: ETF
    R = corr(A);
    val = c'*inv(R)*c;
end