function [R1_t] = dynamic_R1_IR(S_t, T1_0, S_0, TI, eff)
%DYNAMIC_R1_IR Compute dynamic R1 for inversion recovery time-series
%   [R1_t] = dynamic_R1_IR(S_t, T1_0, S_0, TI, eff)
%
% Inputs:
%      S_t - Signal time-series N-d array, with the last dimension time
%
%      T1_0 - Baseline T1, either a scalar or (N-1)-d array, matching the
%      first N-1 dimensions of S_t
%
%      S_0 - Baseline signal, either a scalar or (N-1)-d array, matching the
%      first N-1 dimensions of S_t
%
%      TI - Inversion time, scalar
%
%      eff - Efficiency map for inversion, either a scalar or (N-1)-d array, matching the
%      first N-1 dimensions of S_t. If empty or not set, default 1.0
%
%
% Outputs:
%      R1_t - dynamic R1, N-d array, the same size as S_t
%
%
% Example:
%
% Notes: Non-finite values in R1 (eg due to zero baseline T1) will be set
% to NaN
%
% See also:
%
% Created: 08-Dec-2022
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 5153 
% Copyright: (C) University of Manchester
if ~exist('eff', 'var') || isempty(eff)
    eff = 1;
end

A = 1 - 2*eff .* exp(-TI ./ T1_0);
B = 1 - (S_t ./ S_0).*A;
B(B <= 0) = 0; 

R1_t = -log(0.5*B ./ eff) / TI;

R1_t(~isfinite(R1_t)) = NaN;