clear; clc;
addpath('../../helpers');
load ../../data/spx.mat; 
addpath('../../coint/impl');
addpath('../../coint/deepsearch');
spx_data = csvread('../../data/spx.csv');

twobytwo = 1;
techids  = Extract_Years(pp, twobytwo);
spxids   = [1 508 1016 1524 2032 2540 3048 3556 4064 4572 5080 5588 6096]; %compute with GetPrice(pp_tmp, 236)

for i =flip(1:length(techids)-1)
    start_time  = techids(i); end_time = techids(i+1);
    fprintf('start = %i, end = %i\n', start_time, end_time);
    pp_tmp      = TruncateData(pp, start_time, end_time);
    [ spreads ] = SpreadFinder( pp_tmp, 0.80, 1:1232, 1 );
    fname       = sprintf('spreads_%d_%d.mat', start_time, end_time ); 
    save(fname, 'spreads');
    fname       = sprintf('spx_%d_%d.mat', start_time, end_time );
    spx         = spx_data(spxids(i):spxids(i+1));
    save(fname, 'spx');
end