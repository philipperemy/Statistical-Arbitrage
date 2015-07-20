run('../Init.m');
run('Headers.m');

%ALCOA INC, 236, Materials, AVERY DENNISON CORP, 315, Materials, BANK OF NEW YORK MELLON CORP, 349, Financials, 0.011761, 0.053341, 0.058756, 0.844806, 0.881141, 0.810372, 0.435591, 0.384840, 0.411050

%https://tel.archives-ouvertes.fr/hal-01020405/document p5
%http://fr.mathworks.com/help/econ/unit-root-nonstationarity.html Schwert

% very simple
% determine pmin=0 and pmax=12*(N/100)^0.25
% MonteCarloSteps = N. Calculate N adftest t-stats for each p.
% Average them to have pmax+1 (pmin=0) values.
% Select the lag with the most negative: it indicates the value of lag length that produces the most stationary residuals
% The assumption is that the selected lag length is the optimal value, leading to the best results.

ADFTestOptimalLagFinder(pp, 236, 315, 349);
ADFTestOptimalLagFinder(pp, 713, 949, 1187);

% PP TEST doc very interesting
%     The value of 'model' is determined by the growth characteristics of
%     the time series being tested, and should be chosen with a specific
%     testing strategy in mind. As discussed in Elder & Kennedy [2],
%     including too many regressors results in lost power, while including
%     too few biases the test in favor of the null. In general, if a series
%     is growing, the 'TS' model provides a reasonable trend-stationary
%     alternative to a unit-root process with drift. If a series is not
%     growing, 'AR' and 'ARD' models provide reasonable stationary
%     alternatives to a unit-root process without drift. The 'ARD'
%     alternative has mean c/(1-a); the 'AR' alternative has mean 0.

%This is a nice spread 1000:2000

%Can be a figure
ADFTestOptimalLagFinder(pp, 236, 315, 349, 2000:4500)