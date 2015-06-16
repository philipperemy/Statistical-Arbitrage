function [ indx ] = resampleResidual( w )

    M = length(w);
    indx = zeros(1, M);
    
    % "Repetition counts" (plus the random part, later on):
    Ns = floor(M .* w);

    % The "remainder" or "residual" count:
    R = sum( Ns );

    % The number of particles which will be drawn stocastically:
    M_rdn = M-R;

    % The modified weights:
    Ws = (M .* w - floor(M .* w))/M_rdn;

    % Draw the deterministic part:
    % ---------------------------------------------------
    i=1;
    for j=1:M,
        for k=1:Ns(j),
            indx(i) = j;
            i = i + 1;
        end
    end;

    % And now draw the stocastic (Multinomial) part:
    % ---------------------------------------------------
    Q = cumsum(Ws);
    Q(M)=1; % Just in case...

    while (i<=M),
        sampl = rand;  % (0,1]
        j=1;
        while (Q(j)<sampl),
            j=j+1;
        end;
        indx(i) = j;
        i=i+1;
    end
end
