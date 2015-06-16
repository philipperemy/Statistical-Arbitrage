function [ indx ] = resample( w )

global resampleMethod;
global showResamplePlot;

switch resampleMethod,
    case 0, indx=resampleMultinomial(w);        
    case 1, indx=resampleResidual(w);        
    case 2, indx=resampleStratified(w);        
    case 3, indx=resampleSystematic(w);        
    otherwise,
        error('Invalid value of resampleMethod');
end


if (showResamplePlot),    
	% Draw:
	figure(1);
	subplot(223);
	
	N = length(w);
	Q = cumsum(w);
	
	hold off;
	plot([0;Q],0:N,'r');
	axis([0 1 0 N]);
	
	hold on;
	for i=1:N,
        plot([0 Q(indx(i))],[indx(i) indx(i)],'b');
    end
    title('Inverse CDF');
	
	subplot(224);
	REPETS = hist(indx,1:N);
	barh(REPETS);
	axis([0 max(REPETS) 1 N]);
    title('Repetitions for each sample');
end

