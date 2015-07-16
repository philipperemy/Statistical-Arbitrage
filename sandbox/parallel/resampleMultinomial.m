function [ indx ] = resampleMultinomial( w )
    M = length(w);
    indx = zeros(1, M);
    Q = cumsum(w);
    Q(M)=1; % Just in case...

    i=1;
    while (i<=M),
        sampl = rand;  % (0,1]
        j=1;
        while (Q(j)<sampl),
            j=j+1;
        end;
        indx(i)=j;
        i=i+1;
    end
end