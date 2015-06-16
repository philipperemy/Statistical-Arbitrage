function [] = plotParticlesFromFile( f )

l = load(f);

x(:,:)=l(:,1:3);
w = l(:,4);

g = (w ./ max(w)).^0.5;

for i=1:length(w),
    set(plot(x(i,1),x(i,2)),'Color',[g(i) g(i) g(i)]);        
    hold on;
end

% Fondo negro:
g=get(gcf);
set(g.CurrentAxes,'Color',[0 0 0]);

s = sprintf('mean=(%.02f,%.02f,%.1fdeg), std=(%.02f,%.02f,%.1fdeg) ESS=%f',mean(x(:,1)),mean(x(:,2)),mean(x(:,3))*180/pi,std(x(:,1)),std(x(:,2)),std(x(:,3))*180/pi, ESS(w./sum(w)));
title(s);