function [trends] = CrossingTrends( px1, px2, beg_idx, end_idx )
    spr = px1 - px2;
    trends = zeros(1,1); %Dummy init
    for i = beg_idx:end_idx
        
        if(spr(i-1) > 0 && spr(i) < 0)
            trends(i) = -1;
            continue;
        end
        
        if(spr(i-1) < 0 && spr(i) > 0)
            trends(i) = 1;
            continue;
        end
        
        trends(i) = 0;
    end
end