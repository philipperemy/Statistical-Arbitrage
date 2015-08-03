% efficient. 
% 1 = ids are different
% 0 = at least two ids are alike
function [ ret ] = IdsAreDifferent( id1, id2, id3, id4 )
    ret = ~((id1 == id2) || (id2 == id3) || (id3 == id4) || (id4 == id1) || (id1 == id3) || (id2 == id4)); 
end

