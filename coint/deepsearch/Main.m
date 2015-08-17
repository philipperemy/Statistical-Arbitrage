addpath('../../helpers');
load ../../data/spx.mat; 
addpath('../impl');

i_range         = 235:size(pp.px,2);
[ spreads ]     = SpreadFinder( pp, 0.8, i_range, 1 );