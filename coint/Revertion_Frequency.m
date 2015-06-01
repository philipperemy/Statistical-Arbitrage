function [ crossing_pts, freq_resample ] = Revertion_Frequency( spr, ma_period )
    ma_spr = tsmovavg(spr,'s',ma_period,1);
    
    crossing_pts = CrossingPoints(spr, ma_spr, ma_period+1, length(spr));
    freq = size(crossing_pts, 2);
    dists = abs(diff(crossing_pts));
    
    %http://fr.mathworks.com/help/stats/randsample.htmls
    replacement = true;
    N = 1000;
    dists_distrib = randsample(dists, size(dists, 2)*N, replacement);
    freq_resample = sum(dists_distrib)/mean(dists_distrib)/N;
    
    dur = length(spr) - ma_period - 1;
    freq = freq / dur;
    freq_resample = freq_resample / dur;
end