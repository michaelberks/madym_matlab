function [] = DIBEM_postfit_analysis(subject_dir, subject_idx, visit, varargin)
%DIBEM_POSTFIT_ANALYSIS *Insert a one line summary here*
%   [] = DIBEM_postfit_analysis(varargin)
%
% DIBEM_POSTFIT_ANALYSIS uses the U_PACKARGS interface function
% and so arguments to the function may be passed as name-value pairs
% or as a struct with fields with corresponding names. The field names
% are defined as:
%
% Mandatory Arguments:
%
% Optional Arguments:
%
% Outputs:
%
% Example:
%
% Notes:
%
% See also:
%
% Created: 02-Jul-2020
% Author: Michael Berks 
% Email : michael.berks@manchester.ac.uk 
% Phone : +44 (0)161 275 7669 
% Copyright: (C) University of Manchester 

% Unpack the arguments:
args = u_packargs(varargin, 0, ...
    'subject_id', [],...
    'model_base', 'mdm_analysis_DIBEM_',...
    'dynamic_dir', 'dynamic_cropped',...
    'conc_dir', 'mdm_analysis_DIBEM_f7',...
    'T1_dir', 'mdm_analysis_T1',...
    'model_subscripts', {'f137', 'f37', 'f7'},...
    'model_postscript', [],...
    'visit', 1,...
    'n_vols', 80,...
    'injection_image', 20,...
    'use_linear_base', 1,...
    'model_k', [ 5 6 7],...
    'load_parameters', 0,...
    'load_rois', 0,...
    'do_phys_stats', 0,...
    'do_selection_stats', 0,...
    'do_region_analysis', 0,...
    'do_region_opt_analysis', 0,...
    'roi_mask_paths', [],...
    'roi_names', [],...
    'roi_exclusive_matrix', [],...
    'erode_roi_edges', [],...
    'discard_edge_slices', [],...
    'discard_noisy_voxels', 0,...
    'temporal_noise', [],...
    'roi_analysis', [],...
    'Hct', 0.42,...
    'results_dir', []);
clear varargin;

%% Set-up paths to the data for this subject and some other useful default values
dynamic_dir = [subject_dir args.dynamic_dir '/'];
conc_dir = [subject_dir args.conc_dir '/'];
model_base = [subject_dir args.model_base];

%Create folder for the results
if ~isempty(args.results_dir)
    do_save = true;
    results_dir = [...
        subject_dir args.results_dir '/'];
    create_folder(results_dir);
else
    do_save = false;
end

dibem_models = args.model_subscripts;
n_models = length(dibem_models);

fprintf('\n********************************************************\n');
fprintf('Subject %d, visit %d\n', subject_idx, visit);

%% Load in the DCE concentration time-series

%Load dynamic signals and concentration
[C_t] = get_dyn_vols([conc_dir 'Ct_sig'], args.n_vols);

%Compute time-series error
tse = compute_time_series_error(C_t, args.injection_image+10);
if do_save
    save([results_dir '/' args.subject_id '_tse.mat'],...
        'tse');
end

%**********************************************************************
%%
if args.load_parameters    
    %Load model parameters and convert to physiological form
    [model] = load_nested_DIBEM(...
        model_base, dibem_models, args.model_postscript, 1);

    %Apply AIC analysis
    [aic_selected, sse_selected, aic_weight_maps] = ...
        compute_DIBEM_aic(model, dibem_models, args.model_k, args.n_vols,...
        args.use_linear_base, C_t, args.temporal_noise);
    
    %Get parameters for each selected model
    [model.aic, model.sse] = get_DIBEM_aic_params(...
        model, dibem_models, aic_selected, sse_selected);
    model.aic.weight_maps = aic_weight_maps;
    
    %The model structure isn't changed from this point onwards, so save now
    if do_save
        save([results_dir '/' args.subject_id '_model.mat'],...
            'model');
    end
else
    load([results_dir '/' args.subject_id '_model.mat'],...
        'model');
end
    
%**********************************************************************
%%
roi_analysis = args.roi_analysis;
if args.load_rois
    
    %Load in roi masks
    n_rois = length(args.roi_mask_paths);
    [n_y, n_x, n_z] = size(model.aic.selected);
    roi_masks = false(n_y, n_x, n_z, n_rois);
    
    %Set defaults for modifying the roi
    if length(args.erode_roi_edges) ~= n_rois
        args.erode_roi_edges = zeros(n_rois, 1);
    end
    if length(args.discard_edge_slices) ~= n_rois
        args.discard_edge_slices = zeros(n_rois, 1);
    end
    
    for i_roi = 1:n_rois
        roi_masks(:, :, :, i_roi) = load_img_volume(...
            [subject_dir args.roi_mask_paths{i_roi}]) > 0;
    
        %If selected erode the liver to discount likely noisy boundary voxels
        for i_ero = 1:args.erode_roi_edges(i_roi)
            roi_masks(:, :, :, i_roi) = imerode(...
                roi_masks(:, :, :, i_roi), strel('disk', 1));
        end

        %Discard slices from the top/bottom of the liver 
        if args.discard_edge_slices(i_roi)
            roi_slices = squeeze(any(any(...
                roi_masks(:, :, :, i_roi),1),2));
            slice1 = find(roi_slices,...
                args.discard_edge_slices(i_roi), 'first');
            slice2 = find(roi_slices,...
                args.discard_edge_slices(i_roi), 'last');
            roi_masks(:,:,[slice1 slice2],i_roi) = 0;
        end
    end
    combined_mask = any(roi_masks, 4);
    
    %Reject voxels that are outliers in terms of the temporal noise
    if args.discard_noisy_voxels      
        pct_y = prctile(tse(combined_mask), [25 75]);
        reject_thresh = 1.5*(pct_y(2)-pct_y(1)) + pct_y(2);
        reject_mask = repmat(tse > reject_thresh, 1, 1, 1, n_rois);
        n_rejected = squeeze(...
            sum(sum(sum(reject_mask & roi_masks,3),2),1));
        roi_masks(reject_mask) = 0;
        
        fprintf('%s TSE thresh: %0.2e, rejected:',...
            args.subject_id, reject_thresh);
        fprintf(' %d', n_rejected);
        fprintf('\n\n');
    end
    
    %Finally, apply exclusion criteria
    if ~all(size(args.roi_exclusive_matrix) == n_rois)
        args.roi_exclusive_matrix = false(n_rois, n_rois);
    end
    for i_roi = 1:n_rois
        for j_roi = setdiff(1:n_rois, i_roi)
            if args.roi_exclusive_matrix(i_roi, j_roi)
                roi_i = roi_masks(:,:,:,i_roi);
                roi_i(roi_masks(:,:,:,j_roi)) = 0;
                roi_masks(:,:,:,i_roi) = roi_i;
            end
        end
    end
    
    
    roi_analysis.roi_masks = roi_masks;
else
    n_rois = size(roi_analysis.roi_masks, 4);
end

%**********************************************************************
%%
if args.do_selection_stats
    
    %Compute counts matching each model type
    if args.use_linear_base
        model_names = [dibem_models {'Linear'}];
    else
        model_names = dibem_models;
    end
    n_models_l = length(model_names);
    
    roi_analysis.aic.counts = zeros(n_models_l,n_rois);   
    roi_analysis.sse.counts = zeros(n_models_l,n_rois);
    roi_analysis.n_voxels = zeros(1,n_rois);
    
    for i_roi = 1:n_rois
        roi_mask = roi_masks(:,:,:,i_roi);
        roi_analysis.n_voxels(i_roi) = nnz(roi_mask);
        
        for i_model = 1:n_models_l

            roi_analysis.aic.counts(i_model, i_roi) = sum(...
                model.aic.selected(roi_mask)==i_model);
            roi_analysis.sse.counts(i_model, i_roi) = sum(...
                model.sse.selected(roi_mask)==i_model); 
        end
    end
    
    if length(args.roi_names) ~= n_rois
        args.roi_names = cell(n_rois,1);
    end
    display_DIBEM_model_counts(roi_analysis.aic.counts, model_names, ...
        roi_analysis.n_voxels, args.roi_names, []);   

    if do_save
        save([results_dir '/' args.subject_id '_roi_analysis.mat'],...
            'roi_analysis');
    end
end
%%
if args.do_phys_stats

end
%% ******************************************************************
if args.do_region_analysis
      
    dyn_t = get_dyn_times([dynamic_dir 'dyn_'], args.n_vols);
    
    roi_analysis.aic.C_m = zeros(n_y, n_x, n_z, args.n_vols);
    roi_analysis.sse.C_m = zeros(n_y, n_x, n_z, args.n_vols);
    
    for param_model = 1:n_models
        mod_i = dibem_models{param_model};
        
        aic_preferred_i = model.aic.selected == param_model;
        sse_preferred_i = model.sse.selected(:,:,:,1) == param_model;
        
        analysis_dir = model.(mod_i).analysis_dir;
        C_m = zeros(n_y, n_x, n_z, args.n_vols);
        
        for i_vol = 1:args.n_vols
            C_m_i = load_img_volume(...
                [analysis_dir 'Ct_mod' num2str(i_vol) '.hdr']);
            
            C_m(:,:,:,i_vol) = C_m_i;
            
            C_m_aic_i = roi_analysis.aic.C_m(:,:,:,i_vol);
            C_m_aic_i(aic_preferred_i) = C_m_i(aic_preferred_i);
            roi_analysis.aic.C_m(:,:,:,i_vol) = C_m_aic_i;
            
            C_m_sse_i = roi_analysis.sse.C_m(:,:,:,i_vol);
            C_m_sse_i(sse_preferred_i) = C_m_i(sse_preferred_i);
            roi_analysis.sse.C_m(:,:,:,i_vol) = C_m_sse_i;
        end

        roi_analysis.(mod_i).accumulation_map = ...
            make_DIBEM_accumlation_map(C_m, dyn_t, args.injection_image);
        
        for i_roi = 1:n_rois        
            roi = args.roi_names{i_roi};
            
            %Get results for all voxels
            roi_all = roi_masks(:,:,:,i_roi);
            [roi_vals] = get_roi_vals(...
                model.(mod_i), C_t, C_m, roi_all, args.n_vols);           
            roi_analysis.(mod_i).(roi).all = roi_vals;
            
            %Get results for voxels matched to this model            
            roi_selected = roi_all & (aic_preferred_i);
            [roi_vals] = get_roi_vals(...
                model.(mod_i), C_t, C_m, roi_selected, args.n_vols);           
            roi_analysis.(mod_i).(roi).sel = roi_vals;
            
            %Get results for voxels not-matched to this model            
            roi_notselected = roi_all & (model.aic.selected ~= param_model);
            [roi_vals] = get_roi_vals(...
                model.(mod_i), C_t, C_m, roi_notselected, args.n_vols);           
            roi_analysis.(mod_i).(roi).not = roi_vals;

            
            for mask_model = 1:n_models        
                roi_i = roi_all & (model.aic.selected == mask_model);
                [roi_vals] = get_roi_vals(...
                    model.(mod_i), C_t, C_m, roi_i, args.n_vols);
                roi_analysis.(mod_i).(roi).(dibem_models{mask_model})...
                    = roi_vals;
            end
        end
    end
    
    %Compute accumlation maps for AIC/SSE model concentration
    roi_analysis.aic.accumulation_map = ...
        make_DIBEM_accumlation_map(roi_analysis.aic.C_m, dyn_t, ...
        args.injection_image);
    roi_analysis.sse.accumulation_map = ...
        make_DIBEM_accumlation_map(roi_analysis.sse.C_m, dyn_t, ...
        args.injection_image);
    roi_analysis.signal.accumulation_map = ...
        make_DIBEM_accumlation_map(C_t, dyn_t, ...
        args.injection_image);
    
    roi_analysis.aic.indicator_map = ...
        make_DIBEM_indicator_map(...
        model.aic.active,...
        model.aic.cxm,...
        model.aic.selected, roi_analysis.aic.accumulation_map,...
        roi_analysis.signal.accumulation_map);
    
    roi_analysis.sse.indicator_map = ...
        make_DIBEM_indicator_map(...
        model.sse.active,...
        model.sse.cxm,...
        model.sse.selected(:,:,:,1),...
        roi_analysis.sse.accumulation_map,...
        roi_analysis.signal.accumulation_map);
    
    display_DIBEM_factor_counts(...
        roi_analysis.aic.indicator_map, roi_masks, tse);
    
    if do_save
        save([results_dir '/' args.subject_id '_roi_analysis.mat'],...
            'roi_analysis');
        display_DIBEM_factor_counts(...
            roi_analysis.aic.indicator_map,...
            roi_masks, tse,...
            [results_dir '/' args.subject_id '_indicator.txt']);
    end
end
%%-------------------------------------------------------------------------
% End of main function
%%-------------------------------------------------------------------------
function [roi_vals] = get_roi_vals(model, C_t, C_m, roi, n_vols)

    [dyn_conc] = get_dyn_vals(C_t, n_vols, roi);
    roi_vals.dyn_conc = nanmedian(dyn_conc, 1);

    [model_conc] = get_dyn_vals(C_m, n_vols, roi);
    roi_vals.model_conc = nanmedian(model_conc, 1);

    for i_param = 1:length(model.param_names)
        param_i = model.param_names{i_param};
        param_map = model.(param_i);
        param_vals = param_map(roi);
        param_vals(isnan(param_vals)) = [];

        roi_vals.(param_i) = prctile(param_vals, [0 5 25 50 75 95 100]);
    end
    if isfield(model, 'model_fit')
        param_vals = model.model_fit(roi);
        param_vals(isnan(param_vals)) = [];
        roi_vals.model_fit = prctile(param_vals, [0 5 25 50 75 95 100]);
    end