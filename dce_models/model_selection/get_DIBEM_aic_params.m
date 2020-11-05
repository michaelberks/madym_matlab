function [aic_params, sse_params] = get_DIBEM_aic_params(...
    model, model_names, aic_selected, sse_selected)
%GET_DIBEM_AIC_PARAMS Get parameters for the best fitting model by AIC/SSE
%   [model] = get_DIBEM_aic_params(model, model_names, aic_preferred, sse_preferred)
%
% Inputs:
%      model - model structure returned by COMPUTE_DIBEM_AIC
%
%      model_names - names of the DIBEM model fields
%
%      aic_selected - index array indicating which model selected under
%      AIC criteria
%
%      sse_selected - index array indicating which model selected under
%      SSE criteria
%
%
% Outputs:
%      aic_params - structure containing parameters for each physiological
%      regime as selected by AIC criteria
%
%      sse_params - structure containing parameters for each physiological
%      regime as selected by lowest SSE
%
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 01-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
%

%For each physiological regime, go trhough each nested DIBEM models and get
%the parameters for all voxels where that model is selected under aIC/SSE
%criteria
phys_regimes = {'cxm', 'active'};
opt_param_names{1} = DIBEM_cxm_names();
opt_param_names{2} = DIBEM_active_names();

aic_params.selected = aic_selected;
sse_params.selected = sse_selected;

for i_phys = 1:2
    %Set param names in output structures
    param_names = opt_param_names{i_phys};
    mod_phys = phys_regimes{i_phys};
    aic_params.(mod_phys).param_names = param_names;
    sse_params.(mod_phys).param_names = param_names;
    
    %Loop though each physiological parameter
    for i_param = 1:length(param_names)
        
        %Initialise output params for all voxels from first nested model
        sse_params.(mod_phys).(param_names{i_param}) = ...
            model.(model_names{1}).(mod_phys).(param_names{i_param});
        aic_params.(mod_phys).(param_names{i_param}) = ...
            model.(model_names{1}).(mod_phys).(param_names{i_param});

        %Loop through remaining nested models, swapping voxels where that
        %model is selected
        for i_model = 2:length(model_names)
            mod_i = model_names{i_model};
            
            mask_i = aic_selected(:,:,:,1) == i_model;
            aic_params.(mod_phys).(param_names{i_param})(mask_i) =...
                model.(mod_i).(mod_phys).(param_names{i_param})(mask_i);
            
            mask_i = sse_selected(:,:,:,1) == i_model;
            sse_params.(mod_phys).(param_names{i_param})(mask_i) =...
                model.(mod_i).(mod_phys).(param_names{i_param})(mask_i);
            
        end
    end
end   
