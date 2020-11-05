function [aic_preferred, sse_preferred, aic_weights] = ...
    compute_DIBEM_aic(model, model_names, model_k, n_times,...
    use_linear_base, C_t, temporal_noise)
%COMPUTE_DIBEM_AIC compute AIC weights, and which model is preferred for
%set of nested DIBEM model_names
%   [aic_preferred, sse_preferred, aic_weights] = ...
% compute_DIBEM_aic(model, model_k, use_linear_base)
%
% Inputs:
%      model - model structure returned from LOAD_NESTED_DIBEM
%
%      model_names -
%
%      model_k - number of free model parameters for each model
%
%      n_times - 
%
%      use_linear_base - compute a linear baseline
%
%      C_t -
%
%      temporal_noise -
%
%
% Outputs:
%      aic_preferred - *Insert description of input variable here*
%
%      sse_preferred - *Insert description of input variable here*
%
%      aic_weights - *Insert description of input variable here*
%
%
% Example:
%
% Notes:
%
% See also: LOAD_NESTED_DIBEM
%
% Created: 01-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Compute AIC maps for each model - stack these in a 4-D array so we 
% can easily sort by model
[n_y, n_x, n_z] = size(model.(model_names{1}).model_fit);
n_models = length(model_names);

% Fit a linear model to the data as a baseline in the AIC comparison
if use_linear_base
    aic_maps = inf(n_y, n_x, n_z, n_models+1);
    sse_maps = inf(n_y, n_x, n_z, n_models+1);
    [~,~,aic_maps(:,:,:,n_models+1),sse_maps(:,:,:,n_models+1)] =...  
        fit_baseline_model(C_t, ...
        'input_type', 'conc',...
        'temporal_noise', temporal_noise);
else
    aic_maps = inf(n_y, n_x, n_z, n_models);
    sse_maps = inf(n_y, n_x, n_z, n_models);
end


%Compute AIC for each DIBEM model and stack model fits into container
for i_model = 1:n_models
    model_i = model_names{i_model};
    k = model_k(i_model);
    model.(model_i).k = k;
    N = n_times;
    aic_maps(:,:,:,i_model) = ...
        compute_AIC(model.(model_i).model_fit, N, k, true);
    sse_maps(:,:,:,i_model) = model.(model_i).model_fit;
end 

%Compute preffered model and AIC weights
[aic_preferred, aic_weights] = compute_AIC_weights(aic_maps);
[~, sse_preferred] = sort(sse_maps, 4);