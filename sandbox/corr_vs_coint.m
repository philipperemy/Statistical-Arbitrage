addpath('../coint/impl');

b = cumsum(normrnd(0, 1, 10000, 1));
a = cumsum(normrnd(0, 1, 10000, 1));
Y = [a b];

spr = b - 1.5* a;

N_j = 100;
jci_vec = zeros(1, N_j);
corr_vec = zeros(1, N_j);
adf_vec = zeros(1, N_j);
pp_vec = zeros(1, N_j);
r2_vec = zeros(1, N_j);
for j = 1:N_j
    fprintf('j = %d\n', j);
    tic;
    for i = 1:100
        [h,pValue,stat,cValue,mles] = jcitest(Y,'lags', 0);
    end
    jci_vec(j) = toc;

    tic;
    for i = 1:100
        corr(a,b);
    end
    corr_vec(j) = toc;

    tic;
    for i = 1:100
        [h, p_value_adf] = adftest(spr, 'model','TS','lags', 0);
    end
    adf_vec(j) = toc;

    tic;
    for i = 1:100
        [h,pValue,stat,cValue,reg] = pptest(spr, 'model', 'TS', 'lags', 0);
    end
    pp_vec(j) = toc;

    y = a;
    x1 = b;
    x2 = cumsum(normrnd(0, 1, 10000, 1));

    cor_mat = corr([y x1 x2]);
    c = cor_mat(2:3,1);
    R = corr([x1 x2]);

    tic;
    for i = 1:100
        cor_mat = corr([y x1 x2]);
        c = cor_mat(2:3,1);
        R = corr([x1 x2]);
        h = c'*inv(R)*c;
    end
    r2_vec(j) = toc;
end
