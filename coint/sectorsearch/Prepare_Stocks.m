run('../Init.m');
addpath('../../helpers');
addpath('../impl');

twobytwo            = 1;
year_start_dates    = Extract_Years(pp, twobytwo);
last_start          = year_start_dates(1);
last_end            = year_start_dates(end);

[pp_train, pp_test] = Create_Sets( pp, last_start, last_end, 0.7 );
clear pp; %just to avoid any bugs

[stocks, stocks_tst]=   Prepare_Set( pp_train, pp_test );
sectors             =   unique(stocks.sector);
stocks              =   rmfield(stocks, 'px_clean');
stocks              =   rmfield(stocks, 'w');