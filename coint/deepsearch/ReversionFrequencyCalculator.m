function [ freq ] = ReversionFrequencyCalculator( spr, ma_period )
    ma_spr = tsmovavg(spr,'s',ma_period,1);
    crossing_pts = CrossingPointsCalculator(spr - ma_spr, ma_period+1, length(spr));
    freq = length(crossing_points);
	dur = length(spr) - ma_period - 1;
	freq = sum(abs(crossing_pts)) / dur;
end

%http://fr.mathworks.com/help/stats/randsample.htmls
%replacement = true;
%N = 1000;
%dists_distrib = randsample(dists, size(dists, 2)*N, replacement);
%freq_resample = sum(dists_distrib)/mean(dists_distrib)/N;
%freq = freq / dur;
%freq_resample = freq_resample / dur;