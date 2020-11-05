function [model_params] = load_nested_DIBEM(...
    output_dir_root, model_subscripts, model_postscript, convert_to_phys)
%LOAD_NESTED_DIBEM load parameters for a nested set of DIBEM models fitted
%in Madym
%   [model_params] = load_nested_DIBEM(output_dir_root, 
%                   model_subscripts, model_postscript, convert_to_phys)
%
% Inputs:
%      output_dir_root - Path to the model output dirs, up to the
%      subscript for each nested type
%
%      model_subscripts - subscript defining each nested model type,
%      typically specifying which parameters were fixed
%      default {'f137', f37', 'f7'}
%
%      model_postcript - postscript to each model analysis, default []
%
%      convert_to_phys - flag, if true, converts parameters to
%      physiological forms in the active and exchange regimes
%
%
% Outputs:
%      model_params - structure of loaded parameters. Each nested model is
%      one field in the top-level structure, these models are structure
%      that contain maps of the loaded parameters and model residuals/error
%      codes, and if convert_to_phys is true, conversions of the parameters
%      into their physiological derivations
%
%
% Example: [model_params] = load_nested_DIBEM(...
%    '/mdm_analysis_', {'f137', f37', 'f7'}, '_reg', 1)
%       Will load model paramters from the directories
%       /mdm_analysis_f137_reg, /mdm_analysis_f37_reg, /mdm_analysis_f7_reg
%
%       The output structure will contain fields f137, f37 and f7, each of
%       which is structure containing fields:
%       Fpos, Fneg, Kpos, Kneg, f_a, tau_a, model_fit, error_img,
%       active_params, cxm_params and reject_voxels
%       See CONVERT_DIBEM_PARAMS for details of the parameters in the
%       active_params and cxm_params structures
%
% Notes:
%
% See also: CONVERT_DIBEM_PARAMS
%
% Created: 01-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 
%Load in model_params fit and error maps for each model_params type
if ~exist('model_subscripts', 'var') || isempty(model_subscripts)
    model_subscripts = {};
end

if ~exist('model_postscript', 'var')
    model_postscript = {};
end

if ~exist('convert_to_phys', 'var')
    convert_to_phys = true;
end

%Define these - we need to update Madym to be consistent with naming
madym_names = {'Fpos', 'Fneg', 'Kpos', 'Kneg', 'fa', 'aoffset'};
param_names = {'Fpos', 'Fneg', 'Kpos', 'Kneg', 'f_a', 'tau_a'};

%Loop through each model_params subscipt, loading in the parameters
n_models = length(model_subscripts);
model_params = [];
for i_model = 1:n_models
    mod_i = model_subscripts{i_model};

    %Set up path to model_params analysis folder
    analysis_dir = [output_dir_root mod_i model_postscript '\'];
    model_params.(mod_i).analysis_dir = analysis_dir;

    %Load in model_params residuals and error codes
    model_params.(mod_i).model_fit = load_img_volume(...
        [analysis_dir 'residuals.img']);
    model_params.(mod_i).error_img =...
        load_img_volume([analysis_dir 'error_codes.img']);

    %Load model_params parameters
    model_params.(mod_i).param_names = param_names;

    for i_param = 1:length(param_names)
        model_params.(mod_i).(param_names{i_param}) = ...
            load_img_volume([analysis_dir madym_names{i_param}]);
    end       

    %Convert parameters into active and exchange physiological forms
    if convert_to_phys
        [model_params.(mod_i).active,...
            model_params.(mod_i).cxm, model_params.(mod_i).reject_voxels] = ...
            convert_dibem_params(...
            model_params.(mod_i), model_params.(mod_i).model_fit);

        %Set model_params fit to inf for rejected voxels so this model_params won't be
        %selected for those voxels
        model_params.(mod_i).model_fit(...
            model_params.(mod_i).reject_voxels) = Inf; %| model_params.(mod_i).error_img

        %Copy the shared f_a
        for regime = {'cxm', 'active'}
            model_params.(mod_i).(regime{1}).f_a = model_params.(mod_i).f_a;
            model_params.(mod_i).(regime{1}).tau_a = model_params.(mod_i).tau_a;
        end
    end
end 