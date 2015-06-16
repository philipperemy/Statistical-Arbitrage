function [ COV, MEAN ] = covParticlesXY( x,y,w )
% ---------------------------------------------------------
%  function [ COV ] = covParticlesXY( x,y,w )
%
%  Computes the approximated 2x2 covariance matrix from a 
%    set of particles, given their "x" and "y" coordinates and their
%    weights "w". It also computes the average "mean" value.
%  Jose Luis Blanco Claraco, 26-JUN-2006
% ---------------------------------------------------------

% Assure weights are normalized:
w=w./sum(w);

% The mean values:
MEAN(1) = sum(x.*w);
MEAN(2) = sum(y.*w);

% The covariance:
var_x=0;
var_y=0;
var_xy=0;

var_x = sum( ((x - MEAN(1)).^2).*w );
var_y = sum( ((y - MEAN(2)).^2).*w );
var_xy = sum( ((x - MEAN(1)).*(y - MEAN(2))).*w );

COV(1,1)=var_x;
COV(2,2)=var_y;
COV(1,2)=var_xy;
COV(2,1)=var_xy;

