figure

data1 = normrnd(0.43,1,1000,1);
data2 = normrnd(0.45,1,1000,1);

hist(data1,20)
h1 = findobj(gca,'Type','patch');
set(h1,'FaceColor','g','EdgeColor','w','facealpha',0.75)
hold on
hist(data2,20)
h2 = findobj(gca,'Type','patch');
set(h2,'facealpha',0.8);
ylabel('counts')
xlabel('gene-tree-length/species-tree-length')
legend('no HGT','HGT')
title('p = 0.5');

clear;clc;
addpath('..');
addpath('../../helpers');
c = randi([0 1],200,1)*2-1;
d = cumsum(c);
d = d + 1000;
[mid, uppr, lowr ] = bollinger(d, 20, 1, 1.2);
plot([mid uppr lowr d]); grid;
[ pl, balance_cum, trds ] = TestRoulette(d);
plot(balance_cum);