
%portfolio_cumsum_b
%portfolio_cumsum_z

T = length(portfolio_cumsum_b);
disp_seq = 1:T-1;
%FALSE FIGURE THIS OUT
spx_cumsum = ((10000/spx(166))*spx(166:end));
spx_cumsum = spx_cumsum(disp_seq);

Y = [portfolio_cumsum_b(disp_seq)' portfolio_cumsum_z(disp_seq)' spx_cumsum];
plot(Y);
grid on;