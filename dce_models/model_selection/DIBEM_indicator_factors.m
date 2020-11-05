function [factors] = DIBEM_indicator_factors(indicator_value)
%DIBEM_INDICATOR_FACTORS get list of factors from given DIBEM indicator
%score
%   [factors] = DIBEM_indicator_factors(factor_value)
%
% Inputs:
%      indicator_value - DIBEM indicator_value, an integer <= 2^n_factors
%
%
% Outputs:
%      factors - list of factors in this value as cell array of strings
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 30-Jun-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Define list of factors
factors = {
        'i2'
        'i3'
        'i4'
        'accum'
        'accum_signal'
        'vecs > 1'
        've + vp > 1'
        'kef > ki'
        'active params < 0'
        'CXM params < 0'
        'f_a <= 0'
        'f_a >= 1'};
n_factors = length(factors);
    
%set up index container
factors_idx = false(n_factors,1);

%Get power of 2 factors in value
for i_f = n_factors:-1:1
    factor_score = 2^(i_f-1);

    if indicator_value >= factor_score
        factors_idx(i_f) = 1;
        indicator_value = indicator_value - factor_score;
    end
    if ~indicator_value
        break;
    end
end

%Return list of factors as cell of strings
factors = factors(factors_idx);