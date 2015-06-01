clear;
pepsi = csvread('Pepsi.csv');
coca = csvread('Coca.csv');

%pepsi = pepsi(1000:1500);
%coca = coca(1000:1500);

p = polyfit(coca,pepsi,1);
length = size(coca,1);

plot(1:length,pepsi,'b',1:length,p(1)* coca+p(2),'r');

Y = [coca, pepsi];

[h,pValue,stat,cValue,reg] = egcitest(Y,'test',{'t1','t2'});
pValue

Cointegration(coca, pepsi)

a = reg(2).coeff(1);
b = reg(2).coeff(2);
plot(Y*[1;-b]-a);


length = 250;
sizeC = size(pepsi, 1)-length;
for i = 1:sizeC
    n_pepsi = pepsi(i: i+length);
    n_coca = coca(i: i+length);
    coint(i) = Cointegration(n_coca, n_pepsi);
    fprintf('%i %i %f\n', i, i+length, coint(i));
end