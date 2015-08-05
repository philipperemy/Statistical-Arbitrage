
function [ freq ] = ReversionFrequencyCalculator( spr, ma_period )
    ma_spr = tsmovavg(spr, 's', ma_period, 1);
    crossing_pts = CrossingPointsCalculator(spr - ma_spr, ma_period+1, length(spr));
	duration = length(spr) - ma_period - 1;
	freq = sum(abs(crossing_pts)) / duration;
end