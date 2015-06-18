function [ ] = Assert_Range( val, min_val, max_val )
    error_msg = sprintf('assert range [%f] not in [%f, %f]', val, min_val, max_val);
    assert(min_val < val & val < max_val, error_msg);
end

