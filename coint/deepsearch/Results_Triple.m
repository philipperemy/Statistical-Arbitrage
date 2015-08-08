clear;clc;
load('rep_no_sectors_triples.mat');
load('rep_same_sector_triples.mat');
load('rep_two_sectors_triples.mat');
Y = [rep_same_sector_triples' rep_two_sectors_triples' rep_no_sectors_triples'];
plot(Y);

T = length(Y);


plot(Y(:,1), 'LineWidth', 2);
legend('All from same sector', 'Location','southeast');
ylabel('Cumsum of triples satisfying the criteria');
xlabel('Triples sorted by R^2 DESC');
xlim([0 T+1000]);

% subplot(1,2,2);
% plot(1:length(Y), Y(:,2), 'r', 1:length(Y), Y(:,3), 'b');
% legend('Two from same sector', 'No from same sector', 'Location','southeast');
% xlabel('Triples sorted by R^2 DESC');
% ylabel('Cumsum of triples satisfying the criteria');

W = Y;
for i = 2:length(Y)
    W(i,:) = W(i,:) / sum(Y(i,:));
end

figure;
plot(W(5:end,:), 'LineWidth', 2);
legend('All from same sector', 'Two from same sector', 'None from same sector', 'Location','northeast');
ylabel('Cumsum of triples satisfying the criteria');
xlabel('Triples sorted by R^2 DESC');
xlim([0 T+1000]);
