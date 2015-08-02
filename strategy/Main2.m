addpath('../helpers/');
load ../data/spx.mat; 
spx = csvread('../data/spx.csv');


%one year is not enough. Take more.
year_start_ids = Extract_Years(pp);
year_start_ids = year_start_ids([1 2*(1:1:floor(25/2))+1]); %select period of two years
percentage = 2/3;
spread_length = 10;

nstd_range = 1.0:0.1:3.0;
wts_range = 0:1;
wsize_range = 5:1:20;

%The problem is that when we form the spread. The test asset prices may
%mismatch in prices
for i = 1:length(year_start_ids)-2
    fprintf('year %i-%i\n', i, i+1);
    [pp_train, pp_test] = Create_Sets( pp, year_start_ids(i), year_start_ids(i+2), percentage );
    
    spreads = SpreadFinder(pp_train, 0.8, 236:240, 0);
    
    %Cross validation of the spread
    c = 1;
    clear balances; %hack
    for spread = spreads(1:spread_length)
        try
            [ profits, vols, sharpes, I ] = CrossValidationParallel( pp_train, spread, @SimpleTradingStrategy, nstd_range, wts_range, wsize_range );
            I(:,4)          = sharpes;
            I(:,5)          = profits;
            I(:,6)          = vols;

            I               = I(all(~isnan(I),2),:); % for nan - rows
            I               = sortrows(I, 4);

            wsize           = I(end,3);
            wts             = I(end,2);
            nstd            = I(end,1);
            fprintf('spread %i, best parameters are (%f, %f, %f)\n', c, wsize, wts, nstd);

            %ctor
            spread_tst_px = GetPrice(pp_test, spread.tuple(1))*spread.beta(1) + GetPrice(pp_test, spread.tuple(2))*spread.beta(2) + GetPrice(pp_test, spread.tuple(3))*spread.beta(3);
            spread_tst = struct('px', spread_tst_px, 'tuple', spread.tuple, 'beta', spread.beta); 
            %plot(horzcat(spread.px', spread_tst'))
            [~, balance_cum]= SimpleTradingStrategy( pp_test, spread_tst, 1, length(spread_tst.px), struct('wsize', wsize, 'wts', wts, 'nstd', nstd), 0 );
            balances(c,:)   = balance_cum;
            fprintf('spread %i analysed\n', c);
            c = c + 1;
        catch ex
            disp(getReport(ex));
            %usually because dim mismatch
        end
    end
    
    len = size(balances,1);
    w = repmat(1/len, len, 1);
    portfolio = balances'*w;
    plot(portfolio);
    pl = diff(portfolio);
    sharpe = (mean(pl)/std(pl))*sqrt(252);
    fprintf('sharpe = %f\n', sharpe);
end