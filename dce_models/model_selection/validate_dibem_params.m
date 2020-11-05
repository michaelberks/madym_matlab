function [reject_samples] = validate_dibem_params(model_params, model_fit, irf_form)
%VALIDATE_DIBEM_PARAMS get list of voxels with invalid physiological
%parameters for the different DIBEM model forms
%   [reject_samples] = validate_dibem_params(model, irf_form)
%
% Inputs:
%      model - *Insert description of input variable here*
%
%      irf_form - *Insert description of input variable here*
%
%
% Outputs:
%      reject_samples - *Insert description of input variable here*
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
        ~isfield(model_params, 'Fpos') || ...
        ~isfield(model_params, 'Fpos') || ...
        ~isfield(model_params, 'Fpos')    
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

%Model forms 2 and 3 have specific forms of parameters we know are invalid,
%the full form is a bit more complicated...
switch irf_form
    case 2
        %For volume parameter to be <= 1, we need Kneg <= Fneg. This holds
        %true for either 'active' or exchange interpretation
        reject_samples = reject_samples | ...
            model_params.Fneg >  model_params.Kneg;
    case 3
        %For the mono-exponential form, the active form will always return
        %Ve as infinite. For Vecs <=1 in the active form we need the
        %Fp + ki <= Kneg, where Fp + ki = (Fpos+Fneg)^2 / Fneg
        Fpki = (model_params.Fpos + model_params.Fneg).^2 ./...
                model_params.Fneg;
            
        reject_samples = reject_samples |...
                (Fpki >  model_params.Kneg);
    case 4
        %The maths gets a little complicated for the full form, so easiest
        %to convert to active/exchnage parameters and directly check the
        %volume
        [~, v_ecs] = ...
            active_params_model_to_phys(F_pos, F_neg, K_pos, K_neg, 0);
        
    otherwise
        error('irf_form must be 2, 3 or 4');
end


