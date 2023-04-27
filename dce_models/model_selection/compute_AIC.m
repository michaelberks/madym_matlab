function [AIC] = compute_AIC(SSE, N, k, do_correction)
%COMPUTE_AIC computes Aikake Information criterion (AIC) for model fit
%   [AIC] = compute_AIC(SSE, N, k, do_correction)
%
% Inputs:
%      SSE - sum-of-squared errors in model fit, as all other parameters
%      and scalars, this can be an array of any size (eg the model fit map
%      of an image volume)
%
%      N - number of data points in model fit
%
%      k - number of model parameters (this includes the extra variance
%      parameter of the noise term, eg if your model is y = mx + c + e_i,
%      where e_i is IID Gaussian noise, k = 3 (m, c and sigma(e_i)) )
%
%      do_correction - flag to apply small sample correction, if not set,
%      the default value is to use correction if N/k < 40
%
%
% Outputs:
%      AIC - Aikake information criterion for the model
%
%
% Example:
%
% Notes: Glatting G, Kletting P, Reske S, Hohl K, Ring C. Choosing the
% optimal fit function: comparison of the Akaike information criterion
% and the F-test. Med Phys 2007;34:4285–4292
%
% See also: COMPUTE_AIC_WEIGHTS
%
% Created: 22-Nov-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
if ~exist('do_correction', 'var') || isempty(do_correction)
    do_correction = (N / k) < 40;
end

if do_correction
    correction = 2*k*(k+1) / (N - k - 1);
else
    correction = 0;
end

AIC = real(2*k + N*log(SSE / N) + correction);