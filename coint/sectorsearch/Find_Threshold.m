function [ thres ] = Find_Threshold( dist, number_quads )
    sumX    = sum(dist);
    cumsumX = cumsum(dist);
    backwards = sumX - [0 cumsumX(1:end-1)];
    thres = (1 - sum(flip(backwards) < number_quads)*0.01);
end


%sum(quad_1(Find_Threshold(quad_1, 1000)*100:100))
%sum(quad_2(Find_Threshold(quad_2, 1000)*100:100))
%sum(quad_3(Find_Threshold(quad_3, 1000)*100:100))
%sum(quad_4(Find_Threshold(quad_4, 1000)*100:100))
%sum(quad_5(Find_Threshold(quad_5, 1000)*100:100))
%sum(quad_6(Find_Threshold(quad_6, 1000)*100:100))
%sum(quad_7(Find_Threshold(quad_7, 1000)*100:100))
%sum(quad_8(Find_Threshold(quad_8, 1000)*100:100))
%sum(quad_9(Find_Threshold(quad_9, 1000)*100:100))
%sum(quad_10(Find_Threshold(quad_10, 1000)*100:100))
%sum(quad_11(Find_Threshold(quad_11, 1000)*100:100))