clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');
addpath('./surf/');
addpath('strategy-running/');
addpath('../coint/deepsearch/');
load('../data/spx.mat');

pp = TruncateData(pp, 5000, 5999);
%spreads = SpreadFinder(pp, 0.8, 236:260);
%AMR CORP	238	 Industrials	 CRANE CO	452	 Industrials	 DOVER CORP	498	 Industrials	0.022396	0.182313	0.776812	0.869268	0.911668	0.833231	0.306703	0.490675	0.353043	1.150421
[ h, Spread, ~ ] = SpreadConstructor( [238 452 498], [GetPrice(pp,238) GetPrice(pp,452) GetPrice(pp,498)] );

%Spread = spreads(31);

nstd_range = 0.2:0.1:3.0;
wsize_range = 5:2:120;
wts = 0;

rho1    = 0.999674;
sigma1  = 0.072545;
rho2    = 0.755427;
sigma2  = 0.344350;
beta    = 0.369088;
cor     = -0.853284;
vol_tfsvl = Volatility_TFSVL(Spread, rho1, sigma1, rho2, sigma2, beta, cor);
%vol_tfsvl = 0;

[ profits, vols, sharpes, I ] = CrossValidationParallel( pp, Spread, @SimpleTradingStrategy, nstd_range, wts, wsize_range, vol_tfsvl, 1 );
I(:,4)          = sharpes;
I(:,5)          = profits;
I(:,6)          = vols;

I(isnan(I)) = 0;
A = [I(:,3) I(:,1) I(:,4)];
W = sortrows(A, 3);

wsize           = W(end,1);
nstd            = W(end,2);
[~, balance_cum]= SimpleTradingStrategy( pp, Spread, 1, length(Spread.px), struct('wsize', wsize, 'wts', wts, 'nstd', nstd), vol_tfsvl, 0, 1, @Strategy_Simulator, 10000 );

close all;
figure;
[X,Y,Z] = Draw_Surf(I, nstd_range, wsize_range);
drawnow;
fprintf('balance = %f, sharpe = %f, wsize = %f, nstd = %f\n', balance_cum(end), W(end,3), wsize, nstd);

close all;
ZZ = Plot_Draw_Surf(X,Y,Z);
figure; surf(X,Y,Z);
set(gca, 'CameraPosition', [17.1574 -962.9649   22.5870]);
[px,py] = gradient(ZZ);
figure; quiver(px,py);
quiver(px,py);
set(gca,'XTick', 1:2:length(nstd_range));
set(gca,'XTickLabel', nstd_range(1:2:length(nstd_range)));
set(gca,'YTick', 1:9:length(wsize_range));
set(gca,'YTickLabel', wsize_range(1:9:length(wsize_range)));

[mid, uppr, lowr] = bollinger(Spread.px, 73, 0, 0.9);
Y = [mid uppr lowr Spread.px];
plot(Y);
