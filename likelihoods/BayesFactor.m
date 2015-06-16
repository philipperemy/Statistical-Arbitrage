function [ m ] = BayesFactor( logMarginalLikA, logMarginalLikB )
    log_BF = logMarginalLikA - logMarginalLikB;
    m = 2*log_BF;
end

