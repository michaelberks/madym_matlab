function [data_out] = add_rician_noise(data_in, sigma)
%ADD_RICIAN_NOISE Add Rician noise to a sample of data
%   [data_out] = add_rician_noise(data_in, sigma)
%
% Inputs:
%      data_in - N-d array of data
%
%      sigma - standard deviation of Rician noise
%
%
% Outputs:
%      data_out - N-d array of same size as data_in, with Rician noise
%      added
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 25-Jul-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

data_out = abs(data_in +...
    sigma.*randn(size(data_in)) +...
    sqrt(-1).*sigma.*randn(size(data_in)) );
