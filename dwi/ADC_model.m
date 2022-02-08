function [signals] = ADC_model(B_values, S0, ADC)
%ADC_MODEL *Insert a one line summary here*
%   [signals] = ADC_model(B_values, S0, ADC)
%
% Inputs:
%      B_values - 1-D array of B-values
%
%      S0 - ADC model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%      ADC - ADC model parameter, n-D array or scalar. If array much match
%      numel of other model parameters
%
%
% Outputs:
%      signals - 2D array, numel(B_values) x nSamples, signals computed 
%       from ADC model, 1 sample of B-values per row
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

%Make sure B values are single row and model parameters are single column
%(or scalar)
B_values = B_values(:)';

S0 = S0(:);

ADC = ADC(:);

signals = S0 .* exp(-ADC * B_values);