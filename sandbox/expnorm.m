[f,xi] = ksdensity(randn(10000,1)', -4:0.01:4);


sigma = 1;
R = randn(100000,1);
A1 = exp(sigma/2*R).*R;
[f1,x1] = ksdensity(A1');

sigma = 0.2;
R = randn(100000,1);
A2 = exp(sigma/2*R).*R;
[f2,x2] = ksdensity(A2');

sigma = 0.4;
R = randn(100000,1);
A4 = exp(sigma/2*R).*R;
[f4,x4] = ksdensity(A4');

sigma = 0.6;
R = randn(100000,1);
A6 = exp(sigma/2*R).*R;
[f6,x6] = ksdensity(A6');

sigma = 0.8;
R = randn(100000,1);
A8 = exp(sigma/2*R).*R;
[f8,x8] = ksdensity(A8');

plot(x2, f2, x4, f4, x6, f6, x8, f8, x1, f1);
xlim([-5 10]);
legend('sigma(0.2)', 'sigma(0.4)', 'sigma(0.6)', 'sigma(0.8)', 'sigma(1)');
