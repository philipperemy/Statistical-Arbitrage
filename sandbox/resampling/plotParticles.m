function [ ] = plotParticles( x, w )
% Dibuja un conjunto de particulas de "poses" del robot
%
% function [  ] = plotParticles( x, w )
%

N = size(x,2);  
N0 = N - 3; % The first index of the last "triplet"


g = (w ./ max(w)).^0.2;

for i=1:length(w),
    set(plot(x(i,1:3:end),x(i,2:3:end)),'Color',[g(i) g(i) g(i)]);        
    hold on;
end

% Black background:
g=get(gcf);
set(g.CurrentAxes,'Color',[0 0 0]);

COV = covParticlesXY(x(:,N0+1),x(:,N0+2),w);      

COV(1,1)=COV(1,1)+1e-4;
COV(2,2)=COV(2,2)+1e-4;

M = mean(x(:,N0+(1:2)));
error_ellipse(COV,M,'style','r','conf',0.997);
