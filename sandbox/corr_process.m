x = zeros(1,100);
y = zeros(1,100);

y(1) = randn;
x(1) = randn;

for i = 2:100
   x(i) = 0.8 * y(i-1) + sqrt(1-0.8^2) * randn;
   y(i) = randn;
end

lag_1_y = lagmatrix(y',1);
lag_1_y(isnan(lag_1_y)) = 0;
corr(x', lag_1_y) % they are dependent

lag_1_x = lagmatrix(x',1);
lag_1_x(isnan(lag_1_x)) = 0;
corr(y', lag_1_x) % ok they are independent

lag_2_x = lagmatrix(x',2);
lag_2_x(isnan(lag_2_x)) = 0;
corr(y', lag_2_x) % ok they are independent
