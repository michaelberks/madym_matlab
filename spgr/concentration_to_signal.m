function [S_t, S0_hat, R1_t, e_t, St_hat] = ...
    concentration_to_signal(C_t, FA, TR, T1_0, M0, relax_coeff, use_M0_ratio)
%CONCENTRATION_TO_SIGNAL converts
%   [S_t, S0_hat, R1_t, e_t, St_hat] = ...
%       concentration_to_signal(C_t, FA, TR, T1_0, M0, relax_coeff, use_M0_ratio)
%
% Inputs:
%      C_t (2D-array, nSamples x nTimes) - time-series of dynamic 
%       concentrations, one time-series per column
%
%      FA (scalar) - flip-angle in degrees
%
%      TR (scalar) - repeat time in msecs
%
%      T1_0 (1D array, nTimes x 1) - Baseline T1 in msecs for each time-series
%
%      M_0 (1D array, nTimes x 1) - Either the mean baseline signal to
%      which the output signal will be scaled (using the ratio method) or
%      the M0 value for the T1 fit. One value per time-series
%
%      relax_coeff (scalar - 3.4) - Relaxivity coefficient of tissue in msecs
%
%      use_M0_ratio (integer - 5) - If > 0, specifies the number of
%      time-points over which to average, so that in the output signal,
%      mean(S_t(1:use_M0_ratio)) = M0
%
%
% Outputs:
%      S_t (2D-array, nSamples x nTimes) - dynamic signals
%
%      S0_hat (1D-array, nSamples x 1) - scaling factor applied if 
%           using ratio method
%
%      R1_t (2D-array, nSamples x nTimes) - dynamic R1
%
%      e_t (2D-array, nSamples x nTimes) - exp(-TR .* R1_t)
%
%      St_hat (2D-array, nSamples x nTimes) - unscaled signals
%
%
% Example:
%
% Notes:
%
% See also: SIGNAL_TO_CONCENTRATION
%
% Created: 29-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 

%We fit to columns of time-series
if size(C_t, 2) == 1
    C_t = C_t';
end

% Copyright: (C) University of Manchester 
if ~exist('relax_coeff', 'var') || isempty(relax_coeff)
    relax_coeff = 3.4;
end

%We now set this in msecs to be consistent with Madym
relax_coeff = relax_coeff*1e-3;

if ~exist('use_M0_ratio', 'var') || isempty(use_M0_ratio)
    use_M0_ratio = 5;
end

%Make sure signal inputs are all m x 1 
T1_0 = double(T1_0(:));
M0 = double(M0(:));

[num_voxels,num_time_pts] = size(C_t); %n

%R1 is 1/T1
R1_t = zeros(num_voxels,num_time_pts);
for i_t = 1:num_time_pts
    R1_t(:,i_t) = relax_coeff*C_t(:,i_t) + 1./T1_0;
end

e_t = exp(-TR .* R1_t);
a_t = sind(FA)*(1 - e_t);
b_t = 1 - cosd(FA)*e_t;
St_hat = a_t ./ b_t;

if use_M0_ratio > 0
    %The supplied M0 is actually our target pre-contrast signal S(0), so we can compute
    %M0 as the ratio of S(0) to S(0)_hat
    S0 = M0;
    S0_hat = mean(St_hat(:,1:use_M0_ratio),2);
    M0 = S0 ./ S0_hat;
else
    S0_hat = ones(num_voxels,1);
end

S_t = bsxfun(@times, St_hat, M0);



