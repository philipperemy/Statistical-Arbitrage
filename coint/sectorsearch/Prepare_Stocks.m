run('../Init.m');
addpath('../../helpers');

twobytwo            = 1;
year_start_dates    = Extract_Years(pp, twobytwo);
last_start          = year_start_dates(end-1);
last_end            = year_start_dates(end);

[pp_train_2013, pp_test_2013] = Create_Sets( pp, last_start, last_end, 0.7 );
clear pp; %just to avoid any bugs

stocks              =   Prepare_Set( pp_train_2013 );
sectors             =   unique(stocks.sector);
stocks              =   rmfield(stocks, 'px_clean');
stocks              =   rmfield(stocks, 'w');