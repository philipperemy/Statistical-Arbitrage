load mat.mat;
x1 = mat(:,1);
x2 = mat(:,3);

[X1,X2] = meshgrid(x1,x2);

F = zeros(1,1);
c = 1;
for i = 1:length(X1)
    for j = 1:length(X2)
        F(c) = simple_objective5([X1(i,j), X2(i,j)]);
        c = c + 1;
    end
end
F = reshape(F,length(x2),length(x1));
surf(x1,x2,F);


ObjectiveFunction = @simple_objective5;
X0 = [1.5 0];   % Starting point
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0);