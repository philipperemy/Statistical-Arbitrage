clear;
a = csvread('RDSA.csv');
b = csvread('RDSB.csv');

%500/sec
% tic
%     for i = 1:10000
%         Cointegration(a,b);
%     end
% toc
% 

p = polyfit(a,b,1);
length = size(coca,1);

plot(1:length,pepsi,'b',1:length,p(1)* coca+p(2),'r');

%111/sec
tic
    for i = 1:1000
        Y = [a, b];
        [~,pValue,~,~,~] = egcitest(Y);
    end
toc
disp(pValue);

%50/sec
tic
    for i = 1:1000
        Y = [a, b];
        [~,pValue,~,~,~] = jcitest(Y);
    end
toc
disp(pValue);

Y = [a, b];
[~,pValue,~,~,~] = egcitest(Y);

Y = [a, b];
[~,pValue,~,~,~] = jcitest(Y);