rehash toolboxcache;
x = [1 2];
simple_objective(x);
ObjectiveFunction = @simple_objective4;
X0 = [0.5 0.5];   % Starting point
[x,fval,exitFlag,output] = simulannealbnd(ObjectiveFunction,X0)