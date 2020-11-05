function [preferred_model] = select_dibem_form(model_fit, selection_criteria)
%SELECT_DIBEM_FORM selected preferred model at each voxel depending on
%given statistical criteria and model fit
%   [preferred_model] = select_dibem_form(model_fit, selection_criteria)
%
% Inputs:
%      model_fit - N-d array of model fit residuals, with different models
%      stacked in the N-th dimension
%
%      selection_criteria - {'aic', 'ftest', 'model_fit'} criteria used to
%      select model.
%
%
% Outputs:
%      preferred_model - N-d array, with the model index in sorted order of
%      preference along the last dimension
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
