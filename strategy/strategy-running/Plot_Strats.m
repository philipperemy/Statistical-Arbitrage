
%portfolio_cumsum_b
%portfolio_cumsum_z

LEN = length(portfolio_cumsum_b);
disp_seq = 1:LEN-1;
T = 169;
spx_cumsum = ((10000/spx(T))*spx(T:end));
spx_cumsum = spx_cumsum(disp_seq);

Y = [portfolio_cumsum_b(disp_seq)' portfolio_cumsum_z(disp_seq)' spx_cumsum];
plot(Y);
grid on;

xlabel('Time (days)');
ylabel('Valuation (USD)');