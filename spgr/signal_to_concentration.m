function [C_t, R1_t, e_t, St_hat] = ...
    signal_to_concentration(S_t, FA, TR, T1_0, relax_coeff, num_baseline_pts, m_0)
%CONCENTRATION_TO_SIGNAL converts
%   [S_t, S0_hat, R1_t, e_t, St_hat] = ...
%       concentration_to_signal(C_t, FA, TR, T1_0, M0, relax_coeff, use_M0_ratio)
%
% Inputs:
%      S_t (2D-array, nSamples x nTimes) - time-series of dynamic signals, 
%           one time-series per column
%
%      FA (scalar) - flip-angle in degrees
%
%      TR (scalar) - repeat time in msecs
%
%      T1_0 (1D array, nTimes x 1) - Baseline T1 in msecs for each time-series
%
%      relax_coeff (scalar - 3.4) - Relaxivity coefficient of tissue in msecs
%
%      num_baseline_pts (integer - 5) - If > 0, specifies the number of
%      time-points over which to average, so that in the output signal,
%      mean(C_t(1:num_baseline_pts)) = 0
%
%      m_0 (1D array, nTimes x 1) - Required if not using the ratio method.
%      Scales the concentration according to pre-computed baseline signal
%
%
% Outputs:
%      C_t (2D-array, nSamples x nTimes) - dynamic concentrations
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
% See also: CONCENTRATION_TO_SIGNAL
%
% Created: 29-Mar-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('relax_coeff', 'var') || isempty(relax_coeff)
    relax_coeff = 3.4;
end

%We now set this in msecs to be consistent with Madym
relax_coeff = relax_coeff*1e-3;

if ~exist('num_baseline_pts', 'var') || isempty(num_baseline_pts)
    num_baseline_pts = 5;
end

%We fit to columns of time-series
%if size(S_t, 2) == 1
%    S_t = S_t';
%end

%Make sure signal inputs are all 1 x m 
T1_0 = double(T1_0(:));

if num_baseline_pts
    %Compute scaling factor m_0 from S_0 and T1_0, given the fact C_0 = 0
    if exist('m_0', 'var')
        S_0 = m_0;
    else
        S_0 = mean(S_t(:,1:num_baseline_pts),2);
    end

    e_0 = exp(-TR ./ T1_0);
    a_0 = sind(FA)*(1 - e_0);
    b_0 = 1 - cosd(FA)*e_0;
    m_0 = S_0 .* b_0 ./ a_0;
else
    m_0 = double(m_0(:));
end

St_hat = bsxfun(@rdivide, S_t, m_0);
e_t = (sind(FA) - St_hat) ./ (sind(FA) - St_hat*cosd(FA));
e_t(e_t < 0) = min(e_t(e_t >= 0));

R1_t = -log(e_t) / TR;

C_t = bsxfun(@minus, R1_t, 1 ./ T1_0) / relax_coeff;








