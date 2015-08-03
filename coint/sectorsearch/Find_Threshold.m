function [ thres ] = Find_Threshold( dist, number_quads )
    sumX    = sum(dist);
    cumsumX = cumsum(dist);
    backwards = sumX - [0 cumsumX(1:end-1)];
    thres = (1 - sum(flip(backwards) < number_quads)*0.01);
end

