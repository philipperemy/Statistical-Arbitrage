function [ X, Y, Z ] = Draw_Surf( I, nstd_range, wsize_range )

    I(isnan(I)) = 0;
    A = [I(:,3) I(:,1) I(:,4)];

    [X,Y] = meshgrid(nstd_range, wsize_range);
    Z = zeros(length(wsize_range),length(nstd_range)); %same range
    i = 1;
    for nstd = nstd_range
        j = 1;
       for wsize = wsize_range
           Z(j,i) = Extract_Sharpe(nstd, wsize, A);
           j = j+1;
       end
       i = i+1;
    end

    Plot_Draw_Surf( X, Y, Z );
    
%     %it works
%     surf(X, Y, Z, 'LineStyle', 'none', 'FaceColor', 'interp');
%     colormap(cool);
%     colorbar;
%     camlight right;
%     set(gca, 'CameraPosition', [45 35 9.8]);

    
end

