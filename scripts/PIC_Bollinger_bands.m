load disney;
candle(dis_HIGH(end-20:end), dis_LOW(end-20:end), dis_CLOSE(end-20:end),...
dis_OPEN(end-20:end), 'b');

range = 400:600;
figure;
[mid, uppr, lowr] = bollinger(dis_CLOSE(range));
plot([Fit_NA(mid) Fit_NA(uppr) Fit_NA(lowr)]);
hold on;
candle(dis_HIGH(range), dis_LOW(range), dis_CLOSE(range), dis_OPEN(range), 'black');
legend('mid band', 'upper band', 'lower band', 'stock price');
hold on;
xlim([20 200]);
ylabel('Prices (USD)');
xlabel('Time (Days)');



a = GetPriceArray(pp, [236 237 238 239]);