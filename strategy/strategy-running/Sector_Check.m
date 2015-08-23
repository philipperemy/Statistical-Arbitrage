%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD THE SPREADS %
%%%%%%%%%%%%%%%%%%%%%%

T = length(spreads);
sector_id = zeros(T,1);
load ../../data/spx.mat;
for i = 1:T
    Spread = spreads(i);
    v = pp.sector(Spread.tuple)';
    sector_id(i) = numel(unique(v));
end

clc
mean(sector_id == 1)
mean(sector_id == 2)
mean(sector_id == 3)

sum(sector_id == 1)
sum(sector_id == 2)
sum(sector_id == 3)
