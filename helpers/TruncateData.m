%Keep only rows between start and end
function [ pp_ret ] = TruncateData( pp, dateIdx_start, dateIdx_end )
    pp_ret          = pp; %copy ctor
    pp_ret.dt       = pp.dt(dateIdx_start:dateIdx_end);
    pp_ret.px       = pp.px(dateIdx_start:dateIdx_end, :);
    pp_ret.px_clean = pp.px_clean(dateIdx_start:dateIdx_end, :);
    pp_ret.w        = [];
end

