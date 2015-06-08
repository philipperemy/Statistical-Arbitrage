function [crossing_pts, trends] = CrossingPoints( px1, px2, beg_idx, end_idx )
    c = 1;
    spr = px2 - px1;
    crossing_pts = zeros(1,1); %Dummy init
    trends = zeros(1,1); %Dummy init
    for i = beg_idx:end_idx
        if(spr(i-1) > 0 && spr(i) < 0)
            crossing_pts(c) = i;
            trends(c) = -1;
            c = c + 1;
        else if(spr(i-1) < 0 && spr(i) > 0)
                crossing_pts(c) = i;
                trends(c) = 1;
                c = c + 1;
            end
        end
    end
end