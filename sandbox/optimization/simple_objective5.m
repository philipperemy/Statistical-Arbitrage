function  [sharpe] = simple_objective5(x)
   load mat.mat;
   t = mat;
   vec = find(t(:,1) == x(1));
   rowId = 0;
   for vec_i = vec'
        if(t(vec_i,3) == x(2))
            rowId = vec_i;
            break;
        end
   end
   sharpe = t(rowId, 7);
end
%http://www.math.uri.edu/~bkaskosz/flashmo/graph3d/