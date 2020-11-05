function [active_params, cxm_params, reject_samples] = convert_dibem_params(...
    model_params, model_fit)
%CONVERT_DIBEM_PARAMS convert DIBEM model params to their active and exchange
%derivations, compiling list of invalid voxels
%parameters for the different DIBEM model forms
%   [reject_samples] = validate_dibem_params(model, irf_form)
%
% Inputs:
%      model_params - structure of DIBEM model parameters, must contain
%      fields Fpos, Fneg, Kpos and Kneg each of which should be an N-d
%      array with the same dimensions
%
%      model_fit - residuals of the model fit, N-d array with same
%      dimensions as model parameters
%
%
% Outputs:
%      active_params - structure of parameters for the active derivation,
%      contains fields F_p, v_ecs, k_i, k_ef 
%
%      cxm_params - structure of parameters for the active derivation,
%      contains fields F_p, PS, v_e, v_p
%
%      reject_samples - N-d logical array of indices, 1 where samples have
%      no valid derivation for either physical form
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 29-Jun-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

%Check model_params is valid structure - should have fields Fpos, Fneg,
%Kpos, Kneg
if ~isfield(model_params, 'Fpos') || ...
        ~isfield(model_params, 'Fneg') || ...
        ~isfield(model_params, 'Kpos') || ...
        ~isfield(model_params, 'Kneg')    
	error('model_params structure must contain fiels Fpos, Fneg, Kpos and Kneg');
    
end

%For all model forms, the parameters and model fit should be finite (they
%may return NaN or Inf from madym) and non-negative
reject_samples = ~isfinite(model_fit);
for param = {'Fpos', 'Fneg', 'Kpos', 'Kneg'}
    reject_samples = reject_samples | ...
        ~isfinite(model_params.(param{1})) | ...
        (model_params.(param{1}) < 0);
end

%Convert to active and exchange forms
[active_params.F_p, ...
    active_params.v_ecs,...
    active_params.k_i,...
    active_params.k_ef] = ...
    active_params_model_to_phys(...
        model_params.Fpos, model_params.Fneg, ...
        model_params.Kpos, model_params.Kneg, 0);

[cxm_params.F_p, ...
    cxm_params.PS,...
    cxm_params.v_e,...
    cxm_params.v_p] = ...
    two_cx_params_model_to_phys(...
        model_params.Fpos, model_params.Fneg, ...
        model_params.Kpos, model_params.Kneg, 0);
    
%Check all params are finite and non-negative - do this separately for the
%two derivations, then reject any voxel that is invalid for both. We also
%check volume parameters are <= 1
active_reject = reject_samples;
for param = {'F_p', 'v_ecs', 'k_i', 'k_ef'}
    active_reject = active_reject | ...
        ~isfinite(active_params.(param{1})) | ...
        active_params.(param{1}) < 0;
end
active_reject = active_reject | active_params.v_ecs > 1;

cxm_reject = reject_samples;
for param = {'F_p', 'PS', 'v_e', 'v_p'}
    cxm_reject = cxm_reject | ...
        ~isfinite(cxm_params.(param{1})) | ...
        (cxm_params.(param{1}) < 0);
end
cxm_reject = cxm_reject | (cxm_params.v_e + cxm_params.v_p) > 1;

reject_samples = reject_samples | (active_reject & cxm_reject);
 


