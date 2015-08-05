function [ crossing_points ] = CrossingPointsCalculator( spr, beg_idx, end_idx )
    crossing_points = zeros(1,1); %Dummy init
    for i = beg_idx:end_idx
        
        if(spr(i-1) > 0 && spr(i) < 0)
            crossing_points(i) = -1;
            continue;
        end
        
        if(spr(i-1) < 0 && spr(i) > 0)
            crossing_points(i) = 1;
            continue;
        end
        
        crossing_points(i) = 0;
    end
end