function [ x,w ] = particleFilter( x,w, action, observation )
% A simple SIS with resampling filter
% ------------------------------------------------------------------------------------------------------
% function [ x,w ] = particleFilter( x,w, action, observation )
%
% This function implements one step of the simple SIS with resampling
% particle filter for the especific case of a 2D planar robot and a very simple observation model. 
%  The parameters are:
%   x : Mx(N*3) matrix, with M samples over N-steps in a 3-D path.
%   w : A 1xM or Mx1 vector with the particle importance weights.
%   action: A pose increment vector = (Ax,Ay,Aphi)
%   observation: The actual robot pose (a 3-length vector).
%
% Returns:
%   x : Mx((N+1)*3) matrix, the input with the new pose added to each row.
%   w : 1xM matrix with the new weights.
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

M = size(x,1);

global Beta;
global showResamplePlot;

% Motion model
% ------------------------------------------------
[x] = motionModel( x , action );

% Observation Model
% ------------------------------------------------
Wt = observationLikelihood(x,observation);
w = w .* Wt;


% Normalize the weights:
% ------------------------------------------------
w = w ./ sum(w);

% Resample?
% ------------------------------------------------
ESS_before = ESS(w);
if ( ESS_before < Beta*M ), 
    % Find the indexes for the resample:
    indx = resample(w);
    
    % Duplicate particles & set equal weights:
    x(:,:)=x(indx,:);
    w=ones(M,1)./M;
    
    if (showResamplePlot),
        disp(sprintf('Resampling... ESS: %0.3f -> %0.3f',ESS_before, ESS(w)));
    end
end

