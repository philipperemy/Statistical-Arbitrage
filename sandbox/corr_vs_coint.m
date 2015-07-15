

b = cumsum(normrnd(0, 1, 10000, 1));
a = cumsum(normrnd(0, 1, 10000, 1));
Y = [a b];

spr = b - 1.5* a;

tic;
for i = 1:1000
    [h,pValue,stat,cValue,mles] = jcitest(Y,'lags', 0);
    %c = corr(a,b);
end
toc;

tic;
for i = 1:1000
    c = corr(a,b);
end
toc;

tic;
for i = 1:1000
    [h, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
end
toc;

tic;
for i = 1:1000
    [h,pValue,stat,cValue,reg] = pptest(spr, 'model', 'TS', 'lags', 0);
end
toc;

y = a;
x1 = b;
x2 = cumsum(normrnd(0, 1, 10000, 1));

cor_mat = corr([y x1 x2]);
c = cor_mat(2:3,1);
R = corr([x1 x2]);

c'*inv(R)*c


tic;
for i = 1:1000
    cor_mat = corr([y x1 x2]);
    c = cor_mat(2:3,1);
    R = corr([x1 x2]);
    h = c'*inv(R)*c;
end
toc;