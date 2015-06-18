function [ p_new ] = Parallel_Resample( p, w )
    
    N = length(p);

    %Step 1
    [sort_w, idx] = sort(w); %plot(p(idx), sort_w)
    sort_p = p(idx);
    W = [sort_w sort_p];

    S1 = W(1:2:end,:);
    %k1 = length(S1);

    S2 = W(2:2:end,:);
    %k2 = length(S2);

    c1 = sum(S1(:,1));
    c2 = sum(S2(:,1));

    t = max(c1, c2);
    d = min(c1, c2);
    a = t/d;

    n1 = round(N/(1+a));
    n2 = N - n1;

    %Step 5-1
    Q1 = cumsum(S1(:,1));
    p_new_1 = zeros(1,n1);
    for i = 1:n1

        ui = (c1-0)*rand+0; %runif (0, c1]
        j = 1;

        while (Q1(j) < ui),
            j = j + 1;
        end;

        p_new_1(i) = S1(j,2);
    end

    %Step 5-2
    Q2 = cumsum(S2(:,1));
    p_new_2 = zeros(1,n1);
    for i = 1:n2

        ui = (c2-0)*rand+0; %runif (0, c2]
        j = 1;

        while (Q2(j) < ui),
            j = j + 1;
        end;

        p_new_2(i) = S2(j,2);
    end

    p_new = [p_new_1 p_new_2]';
end

