function [U] = runExample(INTERACTIVE,RESAMPLE_METHOD)
% The main script for running the Particle Filter example
% ------------------------------------------------------------------------------------------------------
% function [U] = runExample(INTERACTIVE,RESAMPLE_METHOD)
%
%   INTERACTIVE =       0/1 : do not/do show graphics and wait for key press
%   RESAMPLE_METHOD =   0: Multinomial Resampling
%                       1: Residual Resampling
%                       2: Stratified Resampling
%                       3: Systematic Resampling
%
% This function runs some steps of a particle filter and return some
% statistics.
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

% --------------------------------------------------------
%   PARAMETERS OF THE SIMULATION
% --------------------------------------------------------
global sensorNoiseStd;
global resampleMethod;
global Beta;
global showResamplePlot;

showResamplePlot=INTERACTIVE;
resampleMethod=RESAMPLE_METHOD;

M = 500;                                    % Number of particles
Beta = 0.5;                                 % Resampling threshold for ESS [0,1]
odometry_bias=[0.00 0.00 0.00];             % Bias of each odometry reading
odometry_std=[0.02 0.02 0.001];             % Standard deviation of the gaussian noise of the x,y,phi components of the odometry
sensorNoiseStd=0.05;                        % Noise for the simulated "global position"-sensor readings
% --------------------------------------------------------
%   END OF "PARAMETERS OF THE SIMULATION"
% --------------------------------------------------------


if (INTERACTIVE),
    close all;
end

% Initial distribution of samples & weights:

% All start at (0,0,0)
x = zeros(M,3);         
% or global localization:
%x = [(rand(M,1)-0.5)*0.2 (rand(M,1)-0.5)*0.2 (rand(M,1)-0.5)*2*pi]; %

w = ones(M,1) ./M;

% Initialize the ground truth robot pose:
x_GT = zeros(1,3);

U(1) = M;    

for t=1:20,
	if (INTERACTIVE),
        % Erase plots:
        figure(1);
        subplot(2,2,1:2);
        hold off;
        % Draw the particles:
        plotParticles(x,w);     % Dibujar.
        title('Samples and ground truth');
	
        % Draw the GT:
        set(plot( x_GT(1:3:end), x_GT(2:3:end), 'b.'),'MarkerSize',30);
        axis equal;
       
        % Wait a key and iterate
        pause;
	end    
    
    
    % Generate the desired (x,y,phi) action for "t":
    action = [1.0 0 0];
    if mod(t,2)==0,
        action(3) = pi/3;
    end
    
    % Move the "real robot":
    x_GT(t*3+1) = x_GT((t-1)*3+1) + action(1).*cos(x_GT((t-1)*3+3)) - action(2).*sin(x_GT((t-1)*3+3));
    x_GT(t*3+2) = x_GT((t-1)*3+2) + action(1).*sin(x_GT((t-1)*3+3)) + action(2).*cos(x_GT((t-1)*3+3));
    x_GT(t*3+3) = x_GT((t-1)*3+3) + action(3);
    
    % Corrupt the action to emulate odometry noise:
    action = action + odometry_bias + odometry_std .* randn(1,3);

    % Observation: The real robot pose plus certain noise:
    observation = [x_GT(t*3+1) x_GT(t*3+2)] + randn(1,2).*sensorNoiseStd;
    
    % Process one step of the particle filter:
    % (resample.m updates subplots 2 & 3)
    [ x,w ] = particleFilter( x,w, action, observation );
            
    % Compute the number of different PATH hypotheses:
    U(t+1) = countPathHypotheses(x);    

    if (INTERACTIVE),
        disp(sprintf('ESS(%2i) = %f',t,ESS(w)));
        disp(sprintf('# of different path hypotheses: %i',U(t+1)));
    end
end
