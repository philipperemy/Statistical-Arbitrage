% WSIZE, NSTD, SHARPE
function [ out ] = Extract_Sharpe( nstd, wsize, A )
    ret = A(A(:,1) == wsize & A(:,2) ==  nstd,:);
    if(size(ret,1)>1)
        ret = ret(1,:);
    end
    out = ret(3);
end

%     try
%         ids1 = find(A(:,1)==wsize);
%         ids2 = find(A(:,2)==nstd);
%         union = intersect(ids1, ids2);
%         sharpe1 = A(union(1),:);
%         if(length(union) == 1)
%             out = sharpe1(3);
%             return;
%         end
%         sharpe2 = A(union(2),:);
%         if(sharpe1(3) > sharpe2(3))
%             out = sharpe1(3);
%         else 
%             out = sharpe2(3);
%         end
%     catch
%         out = 0;
%     end
