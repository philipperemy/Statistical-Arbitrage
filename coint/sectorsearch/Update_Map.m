function [ vec ] = Update_Map( vec, new_R_2 )
    key = floor(round(new_R_2,2)*100);
    if(key == 0) key = 1; end;
    val = vec(key);
    vec(key) = val + 1;
    %fprintf('R2: %f, key %f\n', new_R_2, key);
end
