function [Y_log_SV] = Compute_Log_Likelihood_Generic_From_Filter(varargin)
    if nargin == 1
    end
    filtername = varargin{1};
    varargin = varargin(2:nargin);
    bootstrapFilterHandler = str2func(strcat('TST_Bootstrap_filter_SIR_', filtername));
    [~, ~, ~, p_y] = bootstrapFilterHandler(varargin);
    steps = size(p_y, 2);
    Y_log_SV = Ex2_Estimate_Py( p_y , steps );
    Y_log_SV = Y_log_SV(length(Y_log_SV));
end