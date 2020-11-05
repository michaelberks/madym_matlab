function [sse, model_signals, model_concentrations] = compute_model_sse(...
    dyn_signals, dyn_times, model_params, imaging_params, aif, T1)
%COMPUTE_MODEL_SSE given array of time-series signals and ETM paremeters,
%compute SSE between actual and modelled signal
%   [sse] = compute_model_sse(dyn_signals, model_params, imaging_params, aif, signal_baseline)
%
% Inputs:
%      dyn_signals - Signals of dynamic sequence, (n_t x n_v), where n_t =
%      number of time points, n_v = number of voxels
%
%      dyn_times - times (in s) of each dynmaic nequence (n_t x 1)
%
%      model_params - Struct containing fields
%           Ktrans, Vp, Ve, offset - each should be an n_v length vector
%
%      imaging_params - Struct containing fields
%           FA, TR, relax_coeff, num_precon
%
%      aif - arterial input function (n_t x 1)
%
%      T1 - previously fitted using variable flip-angle (n_v x 1)
%
%
% Outputs:
%      sse - Sum-of-squared error between real and estimated model signals
%       at each voxel
%
%      model_signals - signal at each voxel as estimated by the fitted
%      model
%
%      concentration at each voxel as estimated by the fitted model
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 21-Aug-2017
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Compute model concentrations at each voxel given model paramas
[model_concentrations] = extended_tofts_model(...
    model_params.Ktrans,...
    model_params.Ve,...
    model_params.Vp,...
    model_params.offset, aif, dyn_times);

%Compute baseline signal
signal_baseline =  mean(dyn_signals(:,1:imaging_params.num_precon),2);

%Convert model concentrations to signal given imgaing params, AIF,
%estimated signal baseline and  T1
[model_signals] = concentration_to_signal(model_concentrations,...
    imaging_params.FA, imaging_params.TR, T1, signal_baseline, imaging_params.relax_coeff);

%Compute sum-of- squared error between model signals and original signals
sse = nansum( (dyn_signals - model_signals).^2, 2);
