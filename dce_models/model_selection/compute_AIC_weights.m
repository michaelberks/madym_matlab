function [AIC_preferred, AIC_weights, AIC_min, AIC_delta] = compute_AIC_weights(AIC_vals)
%COMPUTE_AIC_WEIGHTS given a set of AIC values, compute the index of the
%preferred model, and the associated AIC probability weights relative to
%the optimal model
%   [AIC_preferred, AIC_min, AIC_delta, AIC_weights] = compute_AIC_weights(AIC_vals)
%
% Inputs:
%      AIC_vals - set of AIC values, can be an array of any size, in which
%      the multiple models are assumed to be stacked along the last
%      dimension. Eg AICs for N_m models fitted to an Ny x Nx x Nz volume
%      should be passed in as an Ny x Nx x Nz x Nm array
%
%
% Outputs:
%      AIC_preferred - Index of the preferred model, will be an Ny x Nx x Nz
%           array taking values 1,...,Nm
%
%      AIC_weights - probability values associated with each model
%
%      AIC_min - AIC values of the best model
%
%      AIC_delta - delta values for each model
%
% Example:
%
% Notes: Glatting G, Kletting P, Reske S, Hohl K, Ring C. Choosing the
% optimal fit function: comparison of the Akaike information criterion
% and the F-test. Med Phys 2007;34:4285–4292
%
% See also:
%
% Created: 22-Nov-2018
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
model_dim = ndims(AIC_vals);

[AIC_min, AIC_preferred] = min(AIC_vals, [], model_dim);

%Subtract min from AIC for each model
%Could make use of Matlab auto-dimension expansion, but use
%use bsxfun to be compatible with pre-2017 Matlab
AIC_delta = bsxfun(@minus, AIC_vals, AIC_min);

%Compute normalized weights for each model 
AIC_weights = exp(-AIC_delta / 2);

%Make sure we remove NaNs, otherwise these will wipe out good weights in the normalised sum
AIC_weights(isnan(AIC_weights)) = 0; 
AIC_weights = bsxfun(@rdivide, AIC_weights, sum(AIC_weights, model_dim));