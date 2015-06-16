function [ ] = plotParticlesUnweight( x, w )
% Dibuja un conjunto de particulas de "poses" del robot
%
% function [  ] = plotParticles( x, w )
%

plot(x(:,1),x(:,2),'.');


s = sprintf('mean=(%.02f,%.02f,%.1fdeg), std=(%.02f,%.02f,%.1fdeg)',mean(x(:,1)),mean(x(:,2)),mean(x(:,3))*180/pi,std(x(:,1)),std(x(:,2)),std(x(:,3))*180/pi);
title(s);
