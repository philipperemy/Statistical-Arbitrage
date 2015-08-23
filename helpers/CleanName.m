function [ str ] = CleanName( str )
    str = strgsub(char(str), 'UN Equity', '');
    str = strgsub(str, 'UQ Equity', '');
    str = strgsub(str, 'UW Equity', '');
    str = strgsub(str, 'US Equity', '');
end

