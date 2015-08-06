run('Prepare_Stocks.m'); 
load('quad_9.mat');

[ thres ]       = Find_Threshold( quad_9, 10000 ); 
[ quadruples ]  = Dist_Sector_Quad_Retrieve( stocks, sectors(1), 0.25);

quads = quadruples(:,1:4);
quadruples_count = length(quads);
for i = 1:quadruples_count
    quad                = quads(i,:);
    prices              = GetPriceArray( stocks, quad );
    [ h, spread, res ]  = SpreadConstructor4( quad, prices );
    if(h)
       fprintf('spread formed\n');
       prices_tst            = GetPriceArrayTestSet( stocks_tst, quad );
       [ ~, spread_tst, ~ ]  = SpreadConstructor4( quad, prices_tst, spread.beta );
    end
end


%%%scripting for figures
load('quad_1.mat');
load('quad_2.mat');
load('quad_3.mat');
load('quad_4.mat');
load('quad_5.mat');
load('quad_6.mat');
load('quad_7.mat');
load('quad_8.mat');
load('quad_9.mat');
load('quad_10.mat');
load('quad_11.mat');

subplot(2,2,1);
area(quad_2);
subplot(2,2,2);
area(quad_3);
subplot(2,2,3);
area(quad_4);
subplot(2,2,4);
area(quad_4);


subplot(4,3,1);
area(quad_1);
xlabel(sectors(1));
xlim([0 100]);

subplot(4,3,2);
area(quad_2);
xlabel(sectors(2));
xlim([0 100]);

subplot(4,3,3);
area(quad_3);
xlabel(sectors(3));
xlim([0 100]);

subplot(4,3,4);
area(quad_4);
xlabel(sectors(4));
xlim([0 100]);

subplot(4,3,5);
area(quad_5);
xlabel(sectors(5));
xlim([0 100]);

subplot(4,3,6);
area(quad_6);
xlabel(sectors(6));
xlim([0 100]);

subplot(4,3,7);
area(quad_7);
xlabel(sectors(7));
xlim([0 100]);

subplot(4,3,8);
area(quad_8);
xlabel(sectors(8));
xlim([0 100]);

subplot(4,3,9);
area(quad_9);
xlabel(sectors(9));
xlim([0 100]);

subplot(4,3,10);
area(quad_10);
xlabel(sectors(10));
xlim([0 100]);

subplot(4,3,11);
area(quad_11);
xlabel(sectors(11));
xlim([0 100]);
set(gca,'XTick',[1 5 9 10] );
