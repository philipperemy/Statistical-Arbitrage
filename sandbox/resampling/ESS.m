function [ ess ] = ESS( w )
% ------------------------------------------------------------------------------------------------------
% function [ ess ] = ESS( w )
%
% Calcula el "Efective Sample Size". Si es muy bajo es que hay que resamplear las muestras.
% ------------------------------------------------------------------------------------------------------
M = size(w,1);
cv = sum( ( w.*M - 1 ).^2 ) / M;
ess = 1 / ( 1 + cv );
