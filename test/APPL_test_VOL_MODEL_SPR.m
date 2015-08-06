clear; clc;
addpath('../filters/');
addpath('../helpers/');
addpath('../pmcmc/');
addpath('../models/');
addpath('../likelihoods/');
addpath('../coint/deepsearch/');
addpath('../coint/impl/');

%%%% FIRST GET A SPREAD %%%% 09-Sep-2003 - 04-Jun-2006
load('../data/spx.mat');
pp = TruncateData(pp, 5000, 5999);
%spreads = SpreadFinder(pp, 0.8, 236:260);
%AMR CORP	238	 Industrials	 CRANE CO	452	 Industrials	 DOVER CORP	498	 Industrials	0.022396	0.182313	0.776812	0.869268	0.911668	0.833231	0.306703	0.490675	0.353043	1.150421

[ h, spread, ~ ] = SpreadConstructor( [238 452 498], [GetPrice(pp,238) GetPrice(pp,452) GetPrice(pp,498)] );
returns = Compute_Returns(spread.px);
st = SVLogReturns(returns, 0);