function [ indx ] = resampleSystematic( w )

    N = length(w);
    indx = zeros(1, N);
    Q = cumsum(w);

    T = linspace(0,1-1/N,N) + rand/N;
    T(N+1) = 1;

    i=1;
    j=1;

    while (i<=N)
        if (T(i)<Q(j))
            indx(i) = j;
            i = i+1;
        else
            j = j+1;        
       end
    end

end

% function [ indx ] = resampleSystematic( w )
% 
% N = length(w);
% Q = cumsum(w);
% 
% T = linspace(0,1-1/N,N) + rand(1)/N;
% T(N+1) = 1;
% 
% i=1;
% j=1;
% 
% while (i<=N),
%     if (T(i)<Q(j)),
%         indx(i)=j;
%         i=i+1;
%     else
%         j=j+1;        
%     end
% end
% 

