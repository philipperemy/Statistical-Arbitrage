function [ vec ] = GetPrice( pp, i )
    vec = pp.px(:,i);
    vec(isnan(vec)) = [];
end

