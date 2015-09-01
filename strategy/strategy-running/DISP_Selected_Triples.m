function [ ] = DISP_Selected_Triples( pp, spreads, spreads_ids, mat_train, mat_tst )
    tuples = [];
    for i = spreads_ids
       spread = spreads(i);

       c = find(mat_tst(:,5)==i);
       sro = mat_tst(c,3);
       tr = mat_tst(c,4);
       mdd = mat_tst(c,6);
       net = mat_tst(c,1);

       sri = mat_train(i,3);
       s = sprintf('%s & %s & %s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %i \\\\', char(CleanName(pp.tickers(spread.tuple(1)))), char(CleanName(pp.tickers(spread.tuple(2)))), char(CleanName(pp.tickers(spread.tuple(3)))),...
           spread.beta(1), spread.beta(2), spread.beta(3), sri, sro, net, mdd, tr);
       tuples(end+1) = spread.tuple(1);
       tuples(end+1) = spread.tuple(2);
       tuples(end+1) = spread.tuple(3);
       disp(s);
    end

    tuples = unique(tuples);
    for i = 1:length(pp.tickers)
       if(sum(find(i==tuples)~=0))
           s = sprintf('%s & %s & %s \\\\', char(CleanName(pp.tickers(i))), char(pp.names(i)), char(pp.sector(i)));
           disp(s);
       end
    end
end

