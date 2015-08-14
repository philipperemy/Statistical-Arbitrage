function [ RET ] = Plot_Draw_Surf( X, Y, Z, colorBar )

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
    
    RET = Q*coeff_distorsion;
    
    xlabel('Alpha (\alpha)');
    ylabel('Lag (p)');
    zlabel('Sharpe (SR)');
    
    if(exist('colorBar','var'))
        colorbar;
    end;
    colormap(parula(20));
    %     colormap(cool);
%     colorbar;
    %     camlight right;
    set(gca, 'CameraPosition', [16.3457 -372.5140    3.0645]);
    %query is campos
end

