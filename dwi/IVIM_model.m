function [signals] = IVIM_model(B_values, S0, D, f, D_star, full)
%IVIM_MODEL Compute modelled signal
%   [signals] = IVIM_model(B_values, S0, D, f, D_star, full)
%
% Inputs:
%      B_values - 1-D array of B-values
%
%      S0 - IVIM model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%      D - IVIM model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%      f - IVIM model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%      D_star - IVIM model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%      full - bool, default true. Set true if using the ful model. If
%      false, assume D* >> D, so that the 2nd exponential term is ignored
%
%
% Outputs:
%      signals - 2D array, numel(B_values) x nSamples, signals computed 
%       from IVIM model, 1 sample of B-values per row
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Feb-2022
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Assume fitting full model
if ~exist('full', 'var')
    full = 1;
end
if ~full
    %If not using the full model, set D* to infinity, so 2nd term collapses
    %to zero
    D_star = Inf;
end

%Make sure B values are single row and model parameters are single column
%(or scalar)
B_values = B_values(:)';

S0 = S0(:);

D = D(:);

f = f(:);

D_star = D_star(:);

signals = S0 .* (...
    (1 - f) .* exp(-D * B_values) +...
    f .* exp(-D_star * B_values));

