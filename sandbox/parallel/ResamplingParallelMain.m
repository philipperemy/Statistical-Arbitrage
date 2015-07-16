clear;
N = 10000;
mu = 0;
sigma = 1;

p = normrnd(mu, sigma, N, 1);
w = normpdf(p,0,1);
w = w / sum(w);

[sort_p, idx] = sort(p);
plot(sort_p, w(idx)); %this is magic!

%parallel resampling
p1 = p;
for a = 1:1000
    p1 = Parallel_Resample(p1,w);
    fprintf('a = %i\n', a);
end
plot(p1);

p2 = p;
for b = 1:1000
    p2 = p2(resampleMultinomial(w));
    fprintf('b = %i\n', b);
end
plot(p2);