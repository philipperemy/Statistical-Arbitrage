% The main script for obtaining an analysis from PF simulations
% ------------------------------------------------------------------------------------------------------
% The main script for obtaining an analysis from PF simulations
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

close all;
%clear all;

N_EXPERIMENTS=50;

for RESAMPLE_METHOD=0:3,
    for i=1:N_EXPERIMENTS,
        %if mod(i,5)==0, 
            disp(sprintf('METHOD: %i  ->  Iter #%4i/%4i...',RESAMPLE_METHOD,i,N_EXPERIMENTS)); 
        %end;
        U = runExample(0,RESAMPLE_METHOD);
        eval(sprintf('U_%i(i,:)=U;',RESAMPLE_METHOD));
    end
end


figure(1);
subplot(121);

colores(1)='r';
colores(2)='g';
colores(3)='k';
colores(4)='b';

for RESAMPLE_METHOD=0:3,
    eval(sprintf('loglog(mean(U_%i),''%s'');',RESAMPLE_METHOD,colores(RESAMPLE_METHOD+1)));
    hold on;
end

legend('1: Multinomial','2: Residual','3: Stratified','4: Systematic');

xlabel('Timesteps');
title('Number of hypotheses surviving from t=1');
axis([1 size(U_0,2) 1 max(max(U_0))]);


subplot(122);
bar([ mean(U_0(:)) mean(U_1(:)) mean(U_2(:)) mean(U_3(:))]);
xlabel('Method');
title('Expected number of hypotheses surviving from t=1');


w = rand(1000,1);
w = w / sum(w);

ITERATIONS = 10000;
tic;
for i = 1:ITERATIONS
    resampleResidual(w);
end
toc;

%Sample as randsample
% tic;
% for i = 1:ITERATIONS
%     resampleMultinomial(w);
% end
% toc;

%%Simply the best
tic;
for i = 1:ITERATIONS
    resampleStratified(w);
end
toc;

%http://www.site.uottawa.ca/~mbolic/Miodrag_Bolic_files/published/bolic_jasp_2004.pdf
%Said to be the fastest
tic;
for i = 1:ITERATIONS
    resampleSystematic(w);
end
toc;

for i = 1:ITERATIONS
    randsample(1000, 1000, 'true', w);
end
toc;
