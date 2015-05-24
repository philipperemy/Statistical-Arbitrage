function [ vec ] = Convert( pp, i )
    vec = pp.px(:,i);
    vec(isnan(vec)) = [];
end

