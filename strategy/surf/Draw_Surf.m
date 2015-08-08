function [ ] = Draw_Surf( I, nstd_range, wsize_range )

    I(isnan(I)) = 0;
    A = [I(:,3) I(:,1) I(:,4)];

    [X,Y] = meshgrid(nstd_range, wsize_range);
    Z = zeros(96,21); %same range
    i = 1;
    for nstd = nstd_range
        j = 1;
       for wsize = wsize_range
           Z(j,i) = Extract_Sharpe(nstd, wsize, A);
           j = j+1;
       end
       i = i+1;
    end

%     %it works
%     surf(X, Y, Z, 'LineStyle', 'none', 'FaceColor', 'interp');
%     colormap(cool);
%     colorbar;
%     camlight right;
%     set(gca, 'CameraPosition', [45 35 9.8]);

    Q = Z;
    [max_X, max_Y] = size(Z);
    mc = 5;
    for i = 1:mc
        for x = 1:max_X
            for y = 1:max_Y
                xmax = x == max_X;
                ymax = y == max_Y;

                xmin = x == 1;
                ymin = y == 1;

                if(xmax && ymax)
                    Q(x, y) = 0.33*(Q(x-1, y-1)+Q(x-1, y)+Q(x, y-1));
                    continue;
                end

                if(xmin && ymin)
                   Q(x, y) = 0.33*(Q(x+1, y)+Q(x, y+1)+Q(x+1, y+1));
                    continue;
                end

                if(xmax && ymin)
                    Q(x, y) = 0.33*(Q(x, y+1)+Q(x-1, y+1)+Q(x-1, y));
                    continue;
                end

                if(ymax && xmin)
                    Q(x, y) = 0.33*(Q(x+1, y)+Q(x+1, y-1)+Q(x, y-1));
                    continue;
                end

                if(xmax)
                    Q(x, y) = 0.2*(Q(x, y+1)+Q(x-1, y+1)+Q(x-1, y)+ Q(x-1, y-1) + Q(x, y-1));
                    continue;
                end

                if(ymax)
                   Q(x, y) = 0.2*(Q(x-1, y)+Q(x-1, y-1)+Q(x, y-1)+ Q(x+1, y-1) + Q(x+1, y));
                   continue;
                end

                if(xmin)
                    Q(x, y) = 0.2*(Q(x, y+1)+Q(x+1, y+1)+Q(x+1, y)+ Q(x+1, y-1) + Q(x, y-1));
                    continue;
                end

                if(ymin)
                    Q(x, y) = 0.2*(Q(x-1, y)+Q(x-1, y+1)+Q(x, y+1)+ Q(x+1, y+1) + Q(x+1, y));
                    continue;
                end

                Q(x, y) = 0.1250*(Q(x+1, y)+Q(x, y+1)+Q(x-1, y)+Q(x, y-1)+Q(x+1,y+1)+Q(x-1,y-1)+Q(x-1,y+1)+Q(x-1,y+1));
            end
        end
    end

    new_max = max(max(Q));
    real_max = max(max(Z));
    coeff_distorsion = real_max / new_max;
    surf(X,Y,Q*coeff_distorsion);
    colorbar;
    colormap(parula(20));
    %     colormap(cool);
%     colorbar;
    %     camlight right;
    set(gca, 'CameraPosition', [16.3457 -372.5140    3.0645]);
    %query is campos
end

