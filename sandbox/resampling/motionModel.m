function [ X ] = motionModel( X, action )
% ------------------------------------------------------------------------------------------------------
% function [ X ] = Prediction( X, action )
%
% Usa un modelo probabilistico para dar el siguiente estado "x" dado el anterior
%  y la "accion" que se ha realizado. En este caso:
%  x: Es un vector 3x1 de la pose del robot.
%  action: Es un vector 3x1 del cambio de pose que se ha mandado al robot
%  como accion (Ax Ay Aphi)
%
%  J.L. Blanco - University of Malaga, Spain
% ------------------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------------------
% Copyright (c) 2007  Jose Luis Blanco Claraco.
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the "Software"), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% ------------------------------------------------------------------------------------------------------

M = size(X,1);
N = size(X,2);  
N0 = N - 3; % The first index of the last "triplet"

Ax = action(1);
Ay = action(2);
Ap = action(3);

Axy = sqrt(sum(action(1:2).^2));

stdX    = 0.01 * abs(Ax) + 0.01*Axy;
stdY    = 0.01 * abs(Ay) + 0.01*Axy;
stdP    = 0.02 * abs(Ap) + 0.05*Axy;

ccos = cos(X(:,N0+3));
csin = sin(X(:,N0+3));

Ax2 = Ax + randn(M,1) * stdX;
Ay2 = Ay + randn(M,1) * stdY;
Ap2 = Ap + randn(M,1) * stdP;

X(:,N+1) = X(:,N0+1) + Ax2 .* ccos - Ay2 .* csin;
X(:,N+2) = X(:,N0+2) + Ax2 .* csin + Ay2 .* ccos;
X(:,N+3) = X(:,N0+3) + Ap2;
