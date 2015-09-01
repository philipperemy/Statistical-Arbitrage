close all;

data1 = SR_SIMPLE;
data2 = SR_COMPLEX;

for i = 1:(5967-5920)
   SR_COMPLEX(end+1) = 1.1; 
end

SR_COMPLEX = SR_COMPLEX(1:5967);

%data1 = normrnd(0.43,0.1,10000,1);
%data2 = normrnd(0.68,0.1,10000,1);

figure;
hist(data1,20)
h1 = findobj(gca,'Type','patch');
set(h1,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on
hist(data2,20)
h2 = findobj(gca,'Type','patch');
set(h2,'facealpha',0.2);
ylabel('Frequency')
xlabel('Sharpe Ratio')

mu = mean(data1);
plot([mu,mu],ylim,'r--','LineWidth',2);

mu2 = mean(data2);
plot([mu2,mu2],ylim,'b--','LineWidth',2);

legend('Simple','Complex', 'Mean Simple', 'Mean Complex');
xlim([-2.8 4.2]);

figure;
[f, xi] = ksdensity(data1);
plot(xi,f, '--', 'LineWidth', 2, 'Color', 'green');
hold on;
[f, xi] = ksdensity(data2);
plot(xi,f, '--', 'LineWidth', 2, 'Color', 'blue');
legend('SR (Simple Mod)', 'SR (Complex Mod)');
hold off;

figure;
hold on; ecdf(data1); ecdf(data2); hold off;