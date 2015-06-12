% rehash toolboxcache;
% x = [1 2];
% simple_objective(x);
% 
% ObjectiveFunction = @simple_objective;
% X0 = [0.5 0.5];   % Starting point
% [x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0)
% x

x1 = -2:.1:2; x2 = -2:.1:2;
[X1,X2] = meshgrid(x1,x2);

F = zeros(1,1);
c = 1;
for i = 1:length(X1)
    for j = 1:length(X2)
        F(c) = simple_objective3(X1(i,j), X2(i,j));
        c = c + 1;
    end
end
F = reshape(F,length(x2),length(x1));
surf(x1,x2,F);




[X,Y] = meshgrid(-8:.5:8);
Z = simple_objective2(X, Y);
surf(X,Y,Z);
colormap hsv
colorbar;