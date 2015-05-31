function [ freq_resample ] = Revertion_Frequency( spr, ma_period )
    ma_spr = tsmovavg(spr,'s',ma_period,1);
    spr_centered = spr - ma_spr;
    c = 1;
    crossing_pts = zeros(1,1);
    for i = ma_period+1:size(spr_centered, 1)
        if(spr_centered(i-1) > 0 && spr_centered(i) < 0)
            crossing_pts(c) = i;
            c = c + 1;
        else if(spr_centered(i-1) < 0 && spr_centered(i) > 0)
                crossing_pts(c) = i;
                c = c + 1;
            end
        end
    end
    freq = size(crossing_pts, 2);
    dists = abs(diff(crossing_pts));
    
    %http://fr.mathworks.com/help/stats/randsample.htmls
    replacement = true;
    N = 1000;
    dists_distrib = randsample(dists, size(dists, 2)*N, replacement);
    freq_resample = sum(dists_distrib)/mean(dists_distrib)/N;
    
    dur = size(spr_centered, 1) - ma_period - 1;
    freq = freq / dur;
    freq_resample = freq_resample / dur;
end